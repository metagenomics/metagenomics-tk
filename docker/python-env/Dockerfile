FROM python:3.8-slim
ADD requirements.txt /requirements.txt
RUN pip install -r /requirements.txt
RUN apt update && apt install -y procps
