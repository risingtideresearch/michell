FROM caddy:alpine
COPY Caddyfile /etc/caddy/Caddyfile
COPY . /site
RUN rm /site/Caddyfile /site/Dockerfile
