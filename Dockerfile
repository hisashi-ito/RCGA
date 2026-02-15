FROM ruby:2.7-slim

WORKDIR /app
COPY . .

ENTRYPOINT ["ruby", "RCGA.rb"]
