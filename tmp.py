from pathlib import Path
import shutil

SOURCE = Path(
    "outputs/pedal_cell_retro/compact/final_valve_screen/UD/DD"
)

DEST = Path(
    "outputs/pedal_cell_retro/compact/frequency"
)

for label, frequency in (
    ("0p4Hz", 0.4),
    ("0p6Hz", 0.6),
):
    out = DEST / label
    out.mkdir(parents=True, exist_ok=True)

    shutil.copy2(
        SOURCE / "machine.basis.json",
        out / "machine.basis.json",
    )

    text = (SOURCE / "study.toml").read_text()

    old = (
        '"name" = "operation.frequency_hz"\n'
        '"unit" = "Hz"\n'
        '"value" = 0.5'
    )

    new = (
        '"name" = "operation.frequency_hz"\n'
        '"unit" = "Hz"\n'
        f'"value" = {frequency}'
    )

    if old not in text:
        raise RuntimeError("frequency parameter not found")

    text = text.replace(old, new, 1)

    (out / "study.toml").write_text(text)
