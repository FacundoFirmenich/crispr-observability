# Local review candidate. Original wrapper only: see metadata/LICENSE_SCOPE.json.
# Official release tag and manifest digest verified 2026-10-28.
FROM python:3.14.10-slim-bookworm@sha256:5dca57ddc56794d9049dee9ee95a7f071b69d3b018368f2ddc61813649ca2fa0

ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
WORKDIR /opt/mp02
COPY code/ code/
COPY tables/ tables/
COPY raw_replay/ raw_replay/
COPY metadata/DERIVED_INPUT_MANIFEST.json metadata/DERIVED_INPUT_MANIFEST.json
COPY metadata/LICENSE_SCOPE.json metadata/LICENSE_SCOPE.json
COPY THIRD_PARTY_NOTICES.csv README.md CITATION.cff ./
COPY docs/REPRODUCTION_SCOPE.md docs/REPRODUCTION_SCOPE.md
USER 10001:10001
CMD ["python", "-B", "code/run_review_checks.py", "--output-dir", "/tmp/mp02-review"]
