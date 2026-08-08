import requests
import uuid

url = "https://api.yourhospital.org/v1/ingest"

headers = {
    "Content-Type": "application/json",
    "X-Request-Id": str(uuid.uuid4()),
    "X-Idempotency-Key": str(uuid.uuid4()),
}

payload = {
    "patient_key": "pk_7f3a1234",
    "event_ts": "2025-12-30T12:34:56+09:00",
    "payload": {"lab_code": "LAC", "value": 3.2, "unit": "mmol/L"},
}

resp = requests.post(
    url,
    json=payload,
    headers=headers,
    cert=("/path/to/client.crt", "/path/to/client.key"),  # mTLS
    verify="/path/to/ca-chain.crt",  # 서버 인증서 검증
    timeout=2.5
)

print(resp.status_code, resp.text)