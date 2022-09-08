FROM python:3.9

ENV RNA /rna
ENV POETRY_VERSION=1.1.13

WORKDIR $RNA

RUN apt-get update
RUN apt-get upgrade -y

RUN apt-get install -y \
    ca-certificates \
    curl \
    gcc \
    git \
    gzip \
    parallel \
    unzip \
    wget

RUN pip install "poetry==$POETRY_VERSION"

COPY poetry.lock pyproject.toml /rna/

RUN poetry config virtualenvs.create false \
  && poetry install --no-interaction --no-ansi
