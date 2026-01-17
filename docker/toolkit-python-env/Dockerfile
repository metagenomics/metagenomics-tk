FROM python:3.14.0a5-windowsservercore-ltsc2025
ADD requirements.txt /requirements.txt
RUN pip install -r /requirements.txt
RUN apt update && apt install -y procps
