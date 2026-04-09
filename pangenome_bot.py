#!/usr/bin/env python3
"""
Soybean Pangenome Telegram Bot
Runs pangenome_pipeline.py and returns variants to the user.
"""

import os
import sys
import asyncio
import logging
from pathlib import Path
from telegram import Update
from telegram.ext import Application, CommandHandler, MessageHandler, filters, ContextTypes

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(name)s] %(levelname)s: %(message)s")
logger = logging.getLogger("pangenome_bot")

# --- Config ---
TELEGRAM_BOT_TOKEN = os.environ.get("TELEGRAM_BOT_TOKEN", "")
PIPELINE_PATH = Path(os.environ.get(
    "PIPELINE_PATH",
    str(Path(__file__).parent / "pangenome_pipeline.py")
))

# --- Conversation state per user ---
# States: None, "waiting_gfa", "waiting_vcf", "running"
user_state: dict[int, dict] = {}

WELCOME = (
    "🌱 *Soybean Pangenome Bot*\n\n"
    "I can find structural variants in chromosomal regions.\n\n"
    "Send me a message like:\n"
    "_Give me the variants in chromosome 6 between 520000 and 570000_\n\n"
    "Or type /help for more info."
)

HELP = (
    "🌱 *How to use this bot*\n\n"
    "1. Tell me a chromosome and region, e.g.:\n"
    "   `Give me variants in chr6 between 520000 and 570000`\n\n"
    "2. I'll ask for your GFA file path\n"
    "3. Then your VCF file path\n"
    "4. I'll run the pipeline and send you the results\n\n"
    "Type /cancel to start over."
)


async def cmd_start(update: Update, context: ContextTypes.DEFAULT_TYPE):
    await update.message.reply_text(WELCOME, )


async def cmd_help(update: Update, context: ContextTypes.DEFAULT_TYPE):
    await update.message.reply_text(HELP, )


async def cmd_cancel(update: Update, context: ContextTypes.DEFAULT_TYPE):
    uid = update.effective_user.id
    user_state.pop(uid, None)
    await update.message.reply_text("Cancelled. Send a new request whenever you're ready.")


def parse_region(text: str):
    """Try to extract start/end positions from a natural language message."""
    import re
    numbers = re.findall(r'\b(\d{4,})\b', text)
    if len(numbers) >= 2:
        return numbers[0], numbers[1]
    return None, None


async def handle_message(update: Update, context: ContextTypes.DEFAULT_TYPE):
    uid = update.effective_user.id
    text = update.message.text.strip()
    state = user_state.get(uid, {})

    # --- Waiting for GFA path ---
    if state.get("step") == "waiting_gfa":
        gfa_path = text.strip()
        if not Path(gfa_path).exists():
            await update.message.reply_text(
                f"⚠️ File not found: `{gfa_path}`\n\nPlease check the path and try again.",
                
            )
            return
        user_state[uid]["gfa"] = gfa_path
        user_state[uid]["step"] = "waiting_vcf"
        await update.message.reply_text(
            "✅ GFA file found.\n\nNow send me the *full path to your VCF or VCF.gz file*:",
            
        )
        return

    # --- Waiting for VCF path ---
    if state.get("step") == "waiting_vcf":
        vcf_path = text.strip()
        if not Path(vcf_path).exists():
            await update.message.reply_text(
                f"⚠️ File not found: `{vcf_path}`\n\nPlease check the path and try again.",
                
            )
            return
        user_state[uid]["vcf"] = vcf_path
        user_state[uid]["step"] = "running"

        gfa = user_state[uid]["gfa"]
        start = user_state[uid]["start"]
        end = user_state[uid]["end"]
        workdir = Path.home() / f"pangenome_output_{uid}"
        workdir.mkdir(exist_ok=True)

        await update.message.reply_text(
            f"🔬 Running pipeline...\n\n"
            f"GFA: `{gfa}`\n"
            f"VCF: `{vcf_path}`\n"
            f"Region: {start} — {end}\n\n"
            f"This may take a few minutes ⏳",
            
        )

        # Run pipeline
        cmd = [
            sys.executable, str(PIPELINE_PATH),
            "--gfa", gfa,
            "--vcf", vcf_path,
            "--start", start,
            "--end", end,
            "--workdir", str(workdir),
            "--reset"
        ]
        logger.info(f"Running: {' '.join(cmd)}")

        try:
            proc = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )
            stdout, stderr = await asyncio.wait_for(proc.communicate(), timeout=600)
        except asyncio.TimeoutError:
            user_state.pop(uid, None)
            await update.message.reply_text("⏱️ Pipeline timed out after 10 minutes. Try a smaller region.")
            return
        except Exception as e:
            user_state.pop(uid, None)
            await update.message.reply_text(f"❌ Error running pipeline: {e}")
            return

        user_state.pop(uid, None)

        if proc.returncode != 0:
            err = stderr.decode(errors="replace")[-1000:]
            await update.message.reply_text(f"❌ Pipeline failed:\n```\n{err}\n```", )
            return

        # Read output file
        output_file = workdir / "variants_human_readable.txt"
        if not output_file.exists():
            await update.message.reply_text("⚠️ Pipeline finished but no output file found.")
            return

        result = output_file.read_text()
        if not result.strip():
            await update.message.reply_text("✅ Pipeline complete. No variants found in this region.")
            return

        # Send result — split if too long for Telegram (4096 char limit)
        header = f"✅ *Structural variants found in region {start}–{end}*\n\n"
        max_len = 4000 - len(header)

        if len(result) <= max_len:
            await update.message.reply_text(
                header + f"```\n{result}\n```",
                
            )
        else:
            # Send as file
            await update.message.reply_text(header + "_(Result is large — sending as file)_", )
            await update.message.reply_document(
                document=output_file.open("rb"),
                filename="variants_human_readable.txt",
                caption=f"Variants in region {start}–{end}"
            )
        return

    # --- New request ---
    keywords = ["variant", "chromosome", "chr", "between", "region", "position"]
    if any(k in text.lower() for k in keywords):
        start, end = parse_region(text)
        if start and end:
            user_state[uid] = {"step": "waiting_gfa", "start": start, "end": end}
            await update.message.reply_text(
                f"🌱 Got it — region *{start}* to *{end}*\n\n"
                f"Please send me the *full path to your GFA file*:\n"
                f"_(e.g. /Users/yourname/pangenome_work/chr6.gfa)_",
                
            )
        else:
            await update.message.reply_text(
                "I need a start and end position.\n\n"
                "Example: _Give me variants in chromosome 6 between 520000 and 570000_",
                
            )
        return

    # Default
    await update.message.reply_text(
        "I'm the Soybean Pangenome Bot 🌱\n\n"
        "Tell me a chromosomal region to analyze, e.g.:\n"
        "_Give me variants in chr6 between 520000 and 570000_\n\n"
        "Or type /help",
        
    )


def main():
    if not TELEGRAM_BOT_TOKEN:
        print("Error: TELEGRAM_BOT_TOKEN not set.")
        sys.exit(1)
    if not PIPELINE_PATH.exists():
        print(f"Error: pipeline not found at {PIPELINE_PATH}")
        sys.exit(1)

    logger.info(f"Starting Pangenome Bot")
    logger.info(f"Pipeline: {PIPELINE_PATH}")

    app = Application.builder().token(TELEGRAM_BOT_TOKEN).build()
    app.add_handler(CommandHandler("start", cmd_start))
    app.add_handler(CommandHandler("help", cmd_help))
    app.add_handler(CommandHandler("cancel", cmd_cancel))
    app.add_handler(MessageHandler(filters.TEXT & ~filters.COMMAND, handle_message))

    print("🌱 Soybean Pangenome Bot is running. Press Ctrl+C to stop.")
    app.run_polling()


if __name__ == "__main__":
    main()
