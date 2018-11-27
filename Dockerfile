FROM ahmad88me/fcm:latest

RUN mkdir -p /app
WORKDIR /app


COPY scripts /app/scripts
COPY Makefile /app/
COPY *.cpp /app/
COPY *.h /app/
COPY .git /app/.git

CMD ["sh", "scripts/start.sh"]
