FROM caddy:alpine
COPY explainers/Caddyfile /etc/caddy/Caddyfile
COPY explainers/ /site
RUN rm /site/Caddyfile /site/Dockerfile
