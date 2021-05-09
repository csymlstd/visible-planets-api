FROM node:18
WORKDIR /usr/src/app
COPY package.json .
RUN npm install --quiet --production
COPY . .
USER node
