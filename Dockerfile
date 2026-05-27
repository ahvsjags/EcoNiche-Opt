FROM python:3.11-slim

WORKDIR /workspace/econiche-opt
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY . .
RUN pip install --no-cache-dir -e .
ENV PYTHONPATH=/workspace/econiche-opt/src
CMD ["python", "-m", "econiche_opt.cli", "validate-project", "--mode", "demo"]
