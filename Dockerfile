FROM gradle:8.14.3-jdk17

ARG BUILD_HOME=/cellmeasurement

ENV APP_HOME=$BUILD_HOME
WORKDIR $APP_HOME

COPY --chown=gradle:gradle gradle.properties settings.gradle $APP_HOME/
COPY --chown=gradle:gradle gradle $APP_HOME/gradle
COPY --chown=gradle:gradle app $APP_HOME/app

RUN gradle build --no-daemon

RUN echo '#!/bin/sh' > /cellmeasurement.sh && \
    echo 'cd /cellmeasurement' >> /cellmeasurement.sh && \
    echo 'exec gradle run --no-daemon -Dorg.gradle.native=false -Dorg.gradle.configuration-cache=false --no-build-cache --no-problems-report -x compileGroovy -x classes -x processResources "$@"' >> /cellmeasurement.sh && \
    chmod +x /cellmeasurement.sh

ENTRYPOINT ["/cellmeasurement.sh"]
