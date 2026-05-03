const Alexa = require('ask-sdk-core');
const fetch = require('node-fetch');

const API_BASE = 'https://api.spaceleaf.example/v1';

async function getNextSegment(userId, deviceId, previousSegmentId = null) {
  const res = await fetch(`${API_BASE}/news/next`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      user_id: userId,
      device_id: deviceId,
      request_context: {
        previous_segment_id: previousSegmentId,
        timestamp_utc: new Date().toISOString()
      }
    })
  });
  return res.json();
}

async function explainSegment(userId, segmentId) {
  const url = `${API_BASE}/news/explain?user_id=${encodeURIComponent(
    userId
  )}&segment_id=${encodeURIComponent(segmentId)}`;
  const res = await fetch(url);
  return res.json();
}

async function setMode(userId, mode) {
  const res = await fetch(`${API_BASE}/users/mode`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ user_id: userId, mode })
  });
  return res.json();
}

const LaunchRequestHandler = {
  canHandle(handlerInput) {
    return Alexa.getRequestType(handlerInput.requestEnvelope) === 'LaunchRequest';
  },
  async handle(handlerInput) {
    const speakOutput = 'This is Space LEAF Neutral News. Say, play my news, or set my news mode.';
    return handlerInput.responseBuilder.speak(speakOutput).reprompt(speakOutput).getResponse();
  }
};

const PlayNewsIntentHandler = {
  canHandle(handlerInput) {
    const req = handlerInput.requestEnvelope.request;
    return (
      Alexa.getRequestType(handlerInput.requestEnvelope) === 'IntentRequest' &&
      req.intent.name === 'PlayNewsIntent'
    );
  },
  async handle(handlerInput) {
    const { userId } = handlerInput.requestEnvelope.session.user;
    const deviceId = handlerInput.requestEnvelope.context.System.device.deviceId;

    const segment = await getNextSegment(userId, deviceId, null);

    if (segment.action === 'END_OF_FEED') {
      return handlerInput.responseBuilder
        .speak('You are all caught up. There is no more news right now.')
        .getResponse();
    }

    const speakOutput =
      segment.filtered === true
        ? 'Here is your filtered briefing. Some segments were skipped to respect your preferences.'
        : 'Here is your briefing.';

    return handlerInput.responseBuilder
      .speak(speakOutput)
      .addAudioPlayerPlayDirective(
        'REPLACE_ALL',
        segment.audio_url,
        segment.segment_id,
        0,
        null
      )
      .getResponse();
  }
};

const SetModeIntentHandler = {
  canHandle(handlerInput) {
    const req = handlerInput.requestEnvelope.request;
    return (
      Alexa.getRequestType(handlerInput.requestEnvelope) === 'IntentRequest' &&
      req.intent.name === 'SetModeIntent'
    );
  },
  async handle(handlerInput) {
    const { userId } = handlerInput.requestEnvelope.session.user;
    const modeSlot = handlerInput.requestEnvelope.request.intent.slots.mode;
    const modeValue = modeSlot && modeSlot.value;

    // Map spoken value to internal mode
    const modeMap = {
      'neutral news mode': 'NEUTRAL_NEWS_ONLY',
      'neutral mode': 'NEUTRAL_NEWS_ONLY',
      'kid safe mode': 'KID_SAFE',
      'kids mode': 'KID_SAFE',
      'allow all content': 'ALLOW_ALL',
      'no religion': 'NO_RELIGION',
      'no politics': 'NO_POLITICS'
    };

    const internalMode = modeMap[modeValue && modeValue.toLowerCase()];
    if (!internalMode) {
      return handlerInput.responseBuilder
        .speak("I didn't quite catch that mode. You can say neutral news mode or kid safe mode.")
        .reprompt('Which mode would you like?')
        .getResponse();
    }

    await setMode(userId, internalMode);

    const speakOutput = `Your news mode is now set to ${modeValue}.`;
    return handlerInput.responseBuilder.speak(speakOutput).getResponse();
  }
};

const ExplainFilterIntentHandler = {
  canHandle(handlerInput) {
    const req = handlerInput.requestEnvelope.request;
    return (
      Alexa.getRequestType(handlerInput.requestEnvelope) === 'IntentRequest' &&
      req.intent.name === 'ExplainFilterIntent'
    );
  },
  async handle(handlerInput) {
    const { userId } = handlerInput.requestEnvelope.session.user;
    const sessionAttributes = handlerInput.attributesManager.getSessionAttributes();
    const lastFilteredSegmentId = sessionAttributes.lastFilteredSegmentId;

    if (!lastFilteredSegmentId) {
      return handlerInput.responseBuilder
        .speak('I have not skipped any segments recently.')
        .getResponse();
    }

    const explanation = await explainSegment(userId, lastFilteredSegmentId);

    if (!explanation.filtered) {
      return handlerInput.responseBuilder
        .speak('That segment was not filtered for your profile.')
        .getResponse();
    }

    const reasonText =
      explanation.reasons && explanation.reasons.length
        ? explanation.reasons.map(r => r.description).join(' ')
        : 'It did not match your current preferences.';

    return handlerInput.responseBuilder.speak(reasonText).getResponse();
  }
};

const HelpIntentHandler = {
  canHandle(handlerInput) {
    return (
      Alexa.getRequestType(handlerInput.requestEnvelope) === 'IntentRequest' &&
      handlerInput.requestEnvelope.request.intent.name === 'AMAZON.HelpIntent'
    );
  },
  handle(handlerInput) {
    const speakOutput =
      'I can filter your news to remove unsolicited preaching and fear-based messaging. You can say, set my news to neutral mode, or, play my news.';
    return handlerInput.responseBuilder.speak(speakOutput).reprompt(speakOutput).getResponse();
  }
};

const CancelAndStopIntentHandler = {
  canHandle(handlerInput) {
    return (
      Alexa.getRequestType(handlerInput.requestEnvelope) === 'IntentRequest' &&
      (handlerInput.requestEnvelope.request.intent.name === 'AMAZON.CancelIntent' ||
        handlerInput.requestEnvelope.request.intent.name === 'AMAZON.StopIntent')
    );
  },
  handle(handlerInput) {
    return handlerInput.responseBuilder.speak('Goodbye.').getResponse();
  }
};

exports.handler = Alexa.SkillBuilders.custom()
  .addRequestHandlers(
    LaunchRequestHandler,
    PlayNewsIntentHandler,
    SetModeIntentHandler,
    ExplainFilterIntentHandler,
    HelpIntentHandler,
    CancelAndStopIntentHandler
  )
  .lambda();