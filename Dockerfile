FROM gradle:8.14.0-jdk21

ARG BUILD_HOME=/cellmeasurement

ENV APP_HOME=$BUILD_HOME
WORKDIR $APP_HOME

COPY --chown=gradle:gradle gradle.properties settings.gradle $APP_HOME/
COPY --chown=gradle:gradle gradle $APP_HOME/gradle
COPY --chown=gradle:gradle app $APP_HOME/app

RUN gradle build --no-daemon

ENTRYPOINT ["gradle", "run"]
