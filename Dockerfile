FROM python:3.12-slim

WORKDIR /app

COPY requirements.txt ./
COPY app ./app

ENV PORT=8000
EXPOSE 8000

# 依赖随容器启动安装，然后启动服务
CMD ["sh", "-c", "pip install --no-cache-dir -r requirements.txt && uvicorn app.main:app --host 0.0.0.0 --port ${PORT}"]
