FROM ghcr.io/prefix-dev/pixi:0.40.2 AS build

WORKDIR /app
COPY ./pyproject.toml .
COPY ./pixi.lock .
RUN pixi install --locked -e prod
RUN pixi shell-hook -e prod -s bash > /shell-hook
RUN echo "#!/bin/bash" > /app/entrypoint.sh
RUN cat /shell-hook >> /app/entrypoint.sh
RUN echo 'exec "$@"' >> /app/entrypoint.sh

FROM ubuntu:latest AS production
WORKDIR /app
COPY --from=build /app/.pixi/envs/prod /app/.pixi/envs/prod
COPY --from=build --chmod=0755 /app/entrypoint.sh /app/entrypoint.sh
COPY --chmod=0755 ./data /app/data
COPY ./app_code /app/app_code
COPY ./climepi /app/climepi

RUN mkdir /.cache
RUN chmod 777 /.cache
RUN mkdir .chroma
RUN chmod 777 .chroma

ENTRYPOINT [ "/app/entrypoint.sh" ]
CMD ["sh","/app/app_code/run_cluster_app.sh"]