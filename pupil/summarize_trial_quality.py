"""Export the values underlying S3 Table. The participant CSV is restricted."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
from statistics import median

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("retention_csv", type=Path)
parser.add_argument("audit_dir", type=Path)
parser.add_argument("output_dir", type=Path)
args = parser.parse_args()

def read(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))

quality = json.loads((args.audit_dir / "quality_summary.json").read_text())
retention = read(args.retention_csv)
trials = defaultdict(list)
for row in read(args.audit_dir / "raw_signal_quality_by_trial.csv"):
    if row["retained"] == "True":
        trials[(row["pid"], row["visit"])].append(float(row["fixed_0_to_4138_interpolated_fraction"]))

windows = {
    "fixed_0_to_4138": "Common window 0–4138 ms",
    "early_mean_500_to_850": "Early Mean 500–850 ms",
    "early_peak_500_to_2000": "Early Peak 500–2000 ms",
    "late_2128_to_4138": "Late Sustained 2128–4138 ms",
    "response_minus250_to_plus250": "Dilation at Response RT ±250 ms",
    "onset_to_response": "Net Pupil Change 0–RT",
}
summary = []
for key, label in windows.items():
    for visit in ("1", "2"):
        record = quality[visit][key]
        n = quality[visit]["retained_trials_audited"]
        raw = record["raw_low_percent_quantiles"]
        interp = record["interpolated_percent_quantiles"]
        high = record["interpolated_at_least_50pct_trials"]
        summary.append(dict(window=label, visit=visit, retained=n,
            raw_low_median=raw[2], raw_low_q1=raw[1], raw_low_q3=raw[3],
            interpolation_median=interp[2], interpolation_q1=interp[1], interpolation_q3=interp[3],
            at_least_half_interpolated=high, at_least_half_interpolated_percent=100 * high / n))

participants = []
for row in sorted(retention, key=lambda r: (int(r["pid"]), int(r["visit"]))):
    pid, visit = row["pid"], row["visit"]
    values = trials[(pid, visit)]
    record = {key: row[key] for key in ("pid", "visit", "record_available", "recorded",
              "missed_response_excluded", "retained", "congruent_retained", "incongruent_retained",
              "incorrect_responses_retained")}
    if row["record_available"] == "True":
        assert len(values) == int(row["retained"])
        assert len(values) == int(row["congruent_retained"]) + int(row["incongruent_retained"])
    record["median_interpolation_percent"] = 100 * median(values) if values else ""
    participants.append(record)

args.output_dir.mkdir(parents=True, exist_ok=False)
for filename, rows in (("s3_window_diagnostics.csv", summary),
                       ("s3_participant_diagnostics_RESTRICTED.csv", participants)):
    with (args.output_dir / filename).open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
