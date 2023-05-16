
var endTask = [
  {"eng": "Task ended", 
   "fi": "Task ended"},
]

var buttons = [
  {"eng": "Next",
   "fi": "Seuraava"},
  {"eng": "Previous", 
   "fi": "Edellinen"}
]

var recurring = [
  {"eng": "Welcome MMBB battery",
        "fi": "Tervetuloa MMBB-patteriin"},
  {"eng": "Welcome!",
        "fi": "Tervetuloa!"},
  {"eng": "When you are ready to begin, click Next", 
         "fi": "Kun olet valmis aloittamaan, paina Seuraava"},
  {"eng": "End of the task", 
         "fi": "Tehtävä on päättynyt"},
  {"eng": "End of task 2", 
         "fi": "Tehtävä 2 on päättynyt"},
  {"eng": "End of task 3", 
        "fi": "Tehtävä 3 on päättynyt"},
  {"eng": "Loading",      //New
        "fi": "Ladataan"}, //new
  {"eng": "Continue",     //New
        "fi": "Jatka"},        //New
  {"eng": "When you are ready to continue, click Next", 
      "fi": "Kun olet valmis jatkamaan, paina Seuraava"},
  {"eng": "How easy was this task for you?", 
   "fi": "Kuinka helppo tämä tehtävä oli sinulle?"},
  {"eng": "Very easy", 
   "fi": "Todella helppo"}, 
  {"eng": "Very hard", 
   "fi": "Todella haastava"},
  [{"eng": "Welcome to MMBB",
          "fi": "Tervetuloa MMBB-patteriin"},
  {"eng": "For this task you will need the following items:",
    "fi": "Tätä tehtävää varten tarvitset seuraavat esineet:"},
  {"eng": "A flat surface where you can lay down your phone<br><img id='tableSVG' src='../../images/table.svg'></img>",
    "fi": "Tasainen pinta, jolle voit laskea puhelimen<br><img id='tableSVG' src='../../images/table.svg'></img>"},
  {"eng": "A pocket or a purse<br><img id='pocket' src='../../images/pocket.svg'></img>",
   "fi": "Taskun tai käsilaukun<br><img id='pocket' src='../../images/pocket.svg'></img>"},
  {"eng": "Please, use either wired headphones/earplugs or the phone's loudspeaker. Do not use bluetooth hearing-devices.<br> <img id='noBluetooth' src='../../images/noBluetooth.svg'></img>",
   "fi": "Käytä joko langallisia kuulokkeita tai puhelimen kaiutinta. Älä käytä Bluetooth-laitteita."},
  {"eng": "You can take breaks between each task",
   "fi": "Voit pitää tauon jokaisen tehtävän välillä"}],
  [{"eng": "To ensure the correct functioning of the experiment, it is required that you disable the automatic lock of your phone screen. Tap Next to see instructions.",
    "fi": "Jotta koe onnistuu, sinun on poistettava puhelimen näytön automaattinen lukitus. Näet ohjeet napauttamalla Seuraava."},
  {"eng": "On <b>iPhone</b>:",
    "fi": "<b>iPhonin</b> päällä"},
  {"eng": "Settings <br> <img id='instructionLock' src='../../images/instructions/iosInstruction0.png'></img>",
   "fi": "Asetukset <br> <img id='instructionLock' src='../../images/instructions/iosInstruction0.png'></img>"},
  {"eng": "Display & brightness <br> <img id='instructionLock' src='../../images/instructions/iosInstruction1.png'></img>",
   "fi": "Näyttö & kirkkaus <br> <img id='instructionLock' src='../../images/instructions/iosInstruction1.png'></img>"},
  {"eng": "Auto-lock > Never <br> <img id='instructionLock' src='../../images/instructions/iosInstruction2.png'></img>",
    "fi": "Automaattinen lukitus > Ei koskaan <br> <img id='instructionLock' src='../../images/instructions/iosInstruction2.png'></img>"},
  {"eng": "You can re-enable this setting once the experiment has ended.",
   "fi": "Voit ottaa tämän asetuksen uudelleen käyttöön, kun kokeilu on päättynyt."}],
  [{"eng": "To ensure the correct functioning of the experiment, it is required that you disable the automatic lock of your phone screen. Tap Next to see instructions.",
    "fi": "Jotta koe onnistuu, sinun on poistettava puhelimen näytön automaattinen lukitus. Näet ohjeet napauttamalla Seuraava."},
  {"eng": "On <b>Android</b>:",
    "fi": "<b>Androiding</b> päällä"},
  {"eng": "Settings <br> <img id='instructionLock' src='../../images/instructions/androidInstructions0.png'></img>",
    "fi": "Astetukset <br> <img id='instructionLock' src='../../images/instructions/androidInstructions0.png'></img>"},
  {"eng": "Display <br> <img id='instructionLock' src='../../images/instructions/androidInstructions1.png'></img>",
   "fi": "Näyttö & kirkkaus <br> <img id='instructionLock' src='../../images/instructions/androidInstructions1.png'></img>"},
  {"eng": "Screen timeout > 30/15 minutes <br> <img id='instructionLock' src='../../images/instructions/androidInstructions2.png'></img>",
   "fi": "Näytön lukitus > 30/15 minuuttia <br> <img id='instructionLock' src='../../images/instructions/androidInstructions2.png'></img>"},
  {"eng": "You can re-enable this setting once the experiment has ended.",
    "fi": "Voit ottaa tämän asetuksen uudelleen käyttöön, kun kokeilu on päättynyt."},
  ],
  [{"eng": "Not at all",
    "fi": "En lainkaan"},
   {"eng": "Very much",
    "fi": "Erittäin"},
  ]
]
  
var initialInstructions = [
  [{"eng": "Can you hear this song?<br><br>If yes, adjust the volume to a comfortable hearing level and click Continue.<br><br>If you cannot hear the song, try changing to a different computer/phone or to a different internet browser.",
    "fi": "Kuunletko tämän kappaleen?<br>Jos kuulet, säädä äänenvoimakkuutta miellyttävälle kuulotasolle ja napauta Jatka.<br><br>Jos et kuule kappaletta, yritä vaihtaa toiseen tietokoneeseen/puhelimeen tai toiseen internetselaimeen."},
   {"eng": "Yes",
    "fi": "Kyllä"},
   {"eng": "No",
    "fi": "Ei"},
   {"eng": "Is someone helping you to complete this task?",
    "fi": "Auttaako joku sinua tämän tehtävän suorittamisessa?"},
   {"eng": "Can you hear this song?<br><br>If yes, adjust the volume to a comfortable hearing level and click Continue.<br>Make sure that the 'mute' button on the left side of your phone is disabled. <br> If it still doesn't work change to a different phone or to a different internet browser.<br><br>",
    "fi": "Kuunteletko tämän kappaleen?<br><br>Jos kyllä, säädä äänenvoimakkuus miellyttävälle kuulotasolle ja napsauta Seuraava.<br>Jos käytät iPhonea, poista käytöstä puhelimen vasemmalla puolella oleva mykistyspainike. <br> Jos se ei vieläkään toimi, vaihda toiseen puhelimeen tai toiseen Internet-selaimeen.<br><br>"},
   {"eng": "I don't know",
    "fi": "En tiedä"}
  ]
]

var sharedMeasurements = {
  "age": {
    "fi": "Kuinka vanha olet?",
    "eng": "How old are you?"
  },
  "age": {
    "fi": "Kuinka vanha olet?",
    "eng": "How old are you?"
  },
  "gender": {
    "fi": "Mikä on sukupuolesi?",
    "eng": "What is your gender?"
  },
  "genderOption1": {
    "fi": "nainen",
    "eng": "female"
  },
  "genderOption2": {
    "fi": "mies",
    "eng": "male"
  },
  "genderOption3": {
    "fi": "joku muu",
    "eng": "other"
  },
  "genderOption4": {
    "fi": "en halua kertoa",
    "eng": "prefer not to say"
  },
  "languageSpeaksPrimary": {
    "fi": "Mikä on äidinkielesi?",
    "eng": "What is your primary spoken language?"
  },
  "languageSpeaksSecondary": {
    "fi": "Mikä on toissijainen puhuttu kielesi? (jos on)",
    "eng": "What is your secondary spoken language? (If applicable)"
  },
  "height": {
    "fi": "Mikä on pituutesi (cm)?",
    "eng": "What is your height (cm)?"
  },
  "weight": {
    "fi": "Mikä on painosi (kg)?",
    "eng": "What is your weight (Kg)?"
  },
  "yearsPracticeSinging": {
    "fi": "Kuinka monta vuotta olet harrastanut laulamista?",
    "eng": "How many years have you practiced singing?"
  },
  "oftenPracticeSingingActive": {
    "fi": "Kuinka usein harrastit laulamista aktiivisimmassa vaiheessa?",
    "eng": "How often did you practice singing at the most active stage?"
  },
  "oftenPracticeSingingCurrently": {
    "fi": "Kuinka usein harrastat laulamista nykyään?",
    "eng": "How often do you practice singing currently?"
  },
  "yearsPracticeInstrument": {
    "fi": "Kuinka monta vuotta olet harrastanut soittamista?",
    "eng": "How many years have you practiced playing an instrument?"
  },
  "oftenPracticedInstrumentActive": {
    "fi": "Kuinka usein harrastit soittamista aktiivisimmassa vaiheessa?",
    "eng": "How often did you practice playing an instrument at the most active stage?"
  },
  "oftenPracticeInstrumentCurrently": {
    "fi": "Kuinka usein harrastat soittamista nykyään?",
    "eng": "How often do you practice playing an instrument currently?"
  },
  "yearsPracticeDance": {
    "fi": "Kuinka monta vuotta olet harrastanut tanssimista?",
    "eng": "How many years have you practiced dancing?"
  },
  "oftenPracticedDanceActive": {
    "fi": "Kuinka usein harrastit tanssimista aktiivisimmassa vaiheessa?",
    "eng": "How often did you practice dancing at the most active stage? "
  },
  "oftenPracticeDanceCurrently": {
    "fi": "Kuinka usein harrastat tanssimista nykyään?",
    "eng": "How often do you practice dancing currently?"
  },
  "oftenListenMusic": {
    "fi": "Kuinka usein kuuntelet musiikkia aktiivisesti?",
    "eng": "How often do you actively listen to music?"
  },
  "yearsTraining": {
    "fi": "Kuinka monta vuotta olet saanut opetusta peruskoulun ulkopuolella laulamiseen, soittamiseen, ja/tai tanssimiseen?",
    "eng": "How many years of training outside of primary and lower secondary school have you received for singing, playing an instrument, and/or dancing?"
  },
  "howMusical": {
    "fi": "Kuinka musikaaliseksi koet itsesi? 1 (en lainkaan) – 5 (hyvin)",
    "eng": "How musical do you consider yourself? 1 (not at all) – 5 (very)"
  },
  "musicalGenre": {
    "fi": "Mistä musiikkilajeista pidät eniten? (luettele 3)",
    "eng": "Which musical genres do you prefer the most? (list 3)"
  },
  "whichSongArtistsPrefer": {
    "fi": "Mistä lauluista tai artisteista pidät eniten? (luettele 3)",
    "eng": "Which songs or artists do you prefer the most? (list 3)"
  },
  "leastPreferedGenre": {
    "fi": "Mistä musiikkilajeista pidät vähiten? (luettele 3)",
    "eng": "Which musical genres do you prefer the least? (list 3)"
  },
  "leastPreferedArtists": {
    "fi": "Mistä lauluista tai artisteista pidät vähiten? (luettele 3)",
    "eng": "Which songs or artists do you prefer the least? (list 3)"
  },
  "likeMusicEmotion": {
      "eng": "I like to listen to music that contains emotion.",
      "fi": "Pidän sellaisen musiikin kuuntelusta, joka herättää tunteita."
  },
  "getEmotion": {
      "eng": "I get emotional listening to certain pieces of music.",
      "fi": "Liikutun, kun kuuntelen tiettyjä musiikkikappaleita."
  },
  "becomeTearful": {
      "eng": "I can become tearful or cry when I listen to a melody that I like very much.",
      "fi": "Saatan kyynelehtiä tai itkeä kuunnellessani melodiaa, josta pidän todella paljon."
  },
  "feelChill": {
      "eng": "I sometimes feel chills when I hear a melody that I like.",
      "fi": "Toisinaan menen kananlihalle kuunnellessani melodiaa, josta pidän."
  },
  "dontLikeToDance": {
      "eng": "I don’t like to dance, not even with music I like.",
      "fi": "En pidä tanssimisesta edes mielimusiikkini tahtiin."
  },
  "makesMeDance": {
      "eng": "Music often makes me dance.",
      "fi": "Musiikki saa minut usein tanssimaan."
  },
  "humming": {
      "eng": "I can’t help humming or singing along to music that I like.",
      "fi": "En voi olla hyräilemättä tai laulamatta mukana kun kuulen musiikkia, josta pidän."
  },
  "cantStopTapping": {
      "eng": "When I hear a tune I like a lot I can’t help tapping or moving to its beat.",
      "fi": "Kuullessani musiikkikappaleen, josta pidän todella paljon, en malta olla taputtamatta tai liikkumatta sen tahtiin."
  },
  "keepCompany": {
      "eng": "Music keeps me company when I’m alone.",
      "fi": "Musiikki pitää minulle seuraa yksin ollessani."
  },
  "calmsRelaxes": {
      "eng": "Music calms and relaxes me.",
      "fi": "Musiikki rauhoittaa ja rentouttaa minua."
  },
  "chillOut": {
      "eng": "Music helps me chill out.",
      "fi": "Musiikki auttaa minua rentoutumaan."
  },
  "comfortsMe": {
      "eng": "Music comforts me.",
      "fi": "Musiikki lohduttaa minua."
  },
  "hardlyListen": {
      "eng": "In my free time I hardly listen to music.",
      "fi": "En juurikaan kuuntele musiikkia vapaa-aikanani."
  },
  "informMyself": {
      "eng": "I inform myself about music I like.",
      "fi": "Pidän itseni ajan tasalla makuani vastaavasta musiikista."
  },
  "alwaysLooking": {
      "eng": "I’m always looking for new music.",
      "fi": "Etsin aina uutta musiikkia."
  },
  "moneySpend": {
      "eng": "I spend quite a bit of money on music and related items.",
      "fi": "Käytän melko paljon rahaa musiikkiin ja siihen liittyviin asioihin."
  },
  "shareConnection": {
      "eng": "When I share music with someone, I feel a special connection with that person.",
      "fi": "Kun kuuntelen tai esitän musiikkia toisen ihmisen kanssa, tunnen erityistä yhteyttä häneen."
  },
  "bondOtherPeople": {
      "eng": "Music makes me bond with other people.",
      "fi": "Musiikki saa minut kokemaan yhteenkuuluvuutta muiden kanssa."
  },
  "singWithOthers": {
      "eng": "I like to sing or play an instrument with other people.",
      "fi": "Pidän toisten kanssa laulamisesta tai soittamisesta."
  },
  "connectedPerformers": {
      "eng": "At a concert, I feel connected to the performers and the audience.",
      "fi": "Konserteissa koen yhteyden esiintyjien ja yleisön kanssa."
  },
  "backgroundAtmosphere":{
      "eng": "I usually put background music on to make the atmosphere more pleasant.",
      "fi": "Soitan yleensä taustamusiikkia tehdäkseni tunnelmasta miellyttävämmän."
  },
  "busyBackground": {
      "eng": "When I’m busy around the house and no one else is around, I like to have some music on the background.",
      "fi": "Kun olen yksin tehdessäni kotitöitä, pidän taustamusiikin soitosta."
  },
  "afterRough": {
      "eng": "I listen to music to perk up after a rough day. ",
      "fi": "Kuuntelen musiikkia piristyäkseni rankan päivän jälkeen. "
  },
  "exhaustedListen": {
      "eng": "When I’m exhausted, I listen to music to perk up.",
      "fi": "Kun olen uupunut, kuuntelen musiikkia piristyäkseni."
  },
  "magnificentExperience": {
      "eng": "Music has offered me magnificent experiences.",
      "fi": "Musiikki on tarjonnut minulle upeita kokemuksia."
  },
  "feellWholeBody": {
      "eng": "I want to feel the music in my whole body.",
      "fi": "Haluan tuntea musiikin koko kehossani."
  },
  "forgetWorries": {
      "eng": "For me, music is a way to forget about my worries. ",
      "fi": "Musiikki auttaa minua unohtamaan huoleni. "
  },
  "stressfulThoughts": {
      "eng": "When stressful thoughts keep going round and round in my head, I start to listen to music to get them off my mind.",
      "fi": "Kun ahdistavia ajatuksia pyörii päässäni, alan kuuntelemaan musiikkia päästäkseni niistä eroon."
  },
  "reallyAngry": {
      "eng": "When I’m really angry, I feel like listening to some aggressive music.",
      "fi": "Kun olen todella vihainen, minun tekee mieli kuunnella aggressiivista musiikkia."
  },
  "angrySomeone": {
      "eng": "When I’m angry with someone, I listen to music that expresses my anger.",
      "fi": "Kun olen vihainen jollekin, kuuntelen musiikkia joka ilmaisee vihaani."
  },
  "hardExperiences": {
      "eng": "Music has helped me to get through hard experiences.",
      "fi": "Musiikki on auttanut minua selviytymään kovista kokemuksista."
  },
  "distressedClarify": {
      "eng": "When I’m distressed by something, music helps me to clarify my feelings.",
      "fi": "Kun jokin ahdistaa minua, musiikki auttaa selkeyttämään tunteitani."
  },
  "feelsBadComforts": {
      "eng": "When everything feels bad, music understands and comforts me.",
      "fi": "Kun kaikki tuntuu pahalta, musiikki ymmärtää ja lohduttaa minua."
  },
  "feelSadComfort": {
      "eng": "When I’m feeling sad, listening to music comforts me.",
      "fi": "Musiikin kuunteleminen lohduttaa minua kun olen surullinen."
  },
  "ableJudge": {
    "eng": "I am able to judge whether someone is a good singer or not.",
    "fi": "Pystyn arvioimaan onko joku hyvä laulaja vai ei."
  },
  "spotMistakes": {
    "eng": "I find it difficult to spot mistakes in a performance of a song even if I know the tune. ",
    "fi": "Minusta on vaikea havaita virheitä kappaleen esityksessä vaikka tietäisin kappaleen sävelmän."
  },
  "recognizingFamiliar": {
    "eng": "I have trouble recognizing a familiar song when played in a different way or by a different performer.",
    "fi": "Minusta on hankala tunnistaa tuttu kappale kun se soitetaan eri tavalla tai eri esittäjä soittaa sen."
  },
  "canTellOff": {
    "eng": "I can tell when people sing or play out of time with the beat.",
    "fi": "Huomaan kun ihmiset laulavat tai soittavat väärään tahtiin."
  },
  "canTellOutTune": {
    "eng": "I can tell when people sing or play out of tune.",
    "fi": "Huomaan kun ihmiset laulavat tai soittavat epävireessä."
  },
  "noIdeaTune": {
    "eng": "When I sing, I have no idea whether I'm in tune or not.",
    "fi": "Laulaessani en tiedä laulanko epävireessä vai en."
  },
  "usuallyJoin": {
    "eng": "If somebody starts singing a song I don't know, I can usually join in.",
    "fi": "Jos joku alkaa laulamaan minulle tuntematonta kappaletta, voin yleensä liittyä mukaan lauluun."
  },
  "singFromMemory": {
    "eng": "I can sing or play music from memory.",
    "fi": "Osaan laulaa tai soittaa musiikkia ulkomuistista."
  },
  "ableToHitRightNote": {
    "eng": "I am able to hit the right notes when I sing along with a recording.",
    "fi": "Pystyn osumaan oikeisiin nuotteihin laulaessani äänitteen mukana."
  },
  "notAbleHarmony": {
    "eng": "I am not able to sing in harmony when somebody is singing a familiar tune. ",
    "fi": "En pysty laulamaan harmoniassa, kun joku laulaa tuttua säveltä."
  },
  "singingPublic": {
    "eng": "I don’t like singing in public because I’m afraid that I would sing wrong notes.",
    "fi": "En pidä julkisesti laulamisesta, koska pelkään laulavani vääriä nuotteja."
  },
  "singItMyself": {
    "eng": "After hearing a new song two or three times, I can usually sing it by myself.",
    "fi": "Kuultuani uuden kappaleen kahdesti tai kolmesti, pystyn yleensä laulamaan sen itsekseni."
  },
  "easyControlMovement": {
    "eng": "I find it easy to control my movements.",
    "fi": "Minusta on helppo hallita omia liikkeitäni."
  },
  "easyLearnImitate": {
    "eng": "I find it easy to learn or imitate other people's movements.",
    "fi": "Minusta on helppo oppia tai matkia toisten liikkeitä."
  },
  "danceYes": {
    "eng": "If someone asks me to dance, I usually say yes.",
    "fi": "Sanon yleensä kyllä, jos joku pyytää minua tanssimaan."
  },
  "embarrassingDance": {
    "eng": "I find dancing really embarrassing.",
    "fi": "Tanssiminen on minusta todella noloa."
  },
  "whenDanceBetter": {
    "eng": "When I dance, I feel better.",
    "fi": "Minulle tulee parempi olo kun tanssin."
  },
  "feelHaveToDance": {
    "eng": "Sometimes I feel like I just have to dance.",
    "fi": "Joskus koen että minun on vain tanssittava."
  },
}

var tasksIcons = {
  "movement": {
    "eng": "Movement",
    "fi": "Liike"
  },
  "singing": {
    "eng": "Singing",
    "fi": "Laulaminen"
  },
  "survey": {
    "eng": "Survey",
    "fi": "Kysely"
  },
  "emotion": {
    "eng": "Emotion",
    "fi": "Tunteet"
  }
}

var initialPage = {
  "greetings": {
    "eng": "Hello",
    "fi": "Tervetuloa"
  }
}

var emotionTranslations = {
  "openStatement": {
    "eng": "Music and Emotion Survey",
    "fi": "Kysely musiikista ja tunteista"
  },
  "angry": {
    "eng": "Angry",
    "fi": "Viha"
  },
  "happy": {
    "eng": "Happy",
    "fi": "Ilo"
  },
  "sad": {
    "eng": "Sad",
    "fi": "Suru"
  },
  "tender": {
    "eng": "Tender",
    "fi": "Hellyys"
  },
  "instructions1": [
    {
      "eng": "Next, you will listen to clips of music. Please choose the emotion that the music expresses the most.",
      "fi": "Seuraavaksi kuulet musiikkia. Valitse tunnetila, jota musiikki ilmaisee mielestäsi parhaiten."
    },{
      "eng": "You can practice once to get familiar with the task.",
      "fi": "Voit nyt harjoitella, ennen kuin aloitat varsinaisen kyselyn."
  }],
  "instructions2": [{
    "eng": "Thank you for your response! The practice is over, you can start the task by clicking Next",
    "fi": "Kiitos vastauksista! Harjoitus on nyt ohi, ja voit aloittaa tehtävän painamalla Seuraava"
  }],
  "whichEmotion": {
    "eng": "Which of these two emotions does the music express the most?",
    "fi": "Mitä tunnetilaa musiikki ilmaisee parhaiten?"
  },
  "listen": {
    "eng": "Listen...",
    "fi": "Kuuntele…"
  },
  "opening2": {
    "eng": "Music and Emotion Part 2",
    "fi": "Musiikki ja tunteet Osa 2"
  },
  "instructions3": [{
      "eng": "Next, listen to the clips of music and rate what they express in terms of mood, energy, and strength of emotions. You can start the task by clicking \"Next",
      "fi": "Seuraavaksi kuulet musiikkia. Arvioi, mitä eri kappaleet ilmaisevat tunnetilan, energisyyden ja tunteiden voimakkuuden suhteen. Voit aloittaa tehtävän painamalla \"Seuraava"
    },{
      "eng": "Please rate what the music expresses",
      "fi": "Arvioi, mitä musiikki ilmaisee"
    }],
  "mood": {
    "eng": "Mood",
    "fi": "Tunnelita"
  },
  "energy": {
    "eng": "Energy",
    "fi": "Energisyys"
  },
  "howStrongEmotion": {
    "eng": "How strong are the emotions expressed by the music?",
    "fi": "Kuinka voimakkaasti musiikki ilmaisee tunteita?"
  },
  "howMuchLie": {
    "eng": "How much did you like/dislike the music?",
    "fi": "Kuinka paljon pidit/et pitänyt musiikista?"
  },
  "howFamiliar": {
    "eng": "How familiar did the music sound to you?",
    "fi": "Kuinka tutulta musiikki kuulosti?"
  },
  "veryNegative": {
    "eng": "Very negative",
    "fi": "Todella negatiivinen"
  },
  "neutral": {
    "eng": "Neutral",
    "fi": "Neutraali"
  },
  "veryPositive": {
    "eng": "Very positive",
    "fi": "Todella positiivinen"
  },
  "veryLowEnergy": {
    "eng": "Very low energy",
    "fi": "Ei lankaan energinen"
  },
  "veryHighEnergy": {
    "eng": "Very high energy",
    "fi": "Todella vahvasti energinen"
  },
  "notAtAll": {
    "eng": "Not at all",
    "fi": "Ei lainkaan"
  },
  "veryStrong": {
    "eng": "Very strong",
    "fi": "Todella voimakkaita"
  },
  "stronglyDisliked": {
    "eng": "Strongly disliked",
    "fi": "En pidä lainkaan"
  },
  "neutral": {
    "eng": "Neutral",
    "fi": "Neutraali"
  },
  "choose": {
    "eng": "Choose",
    "fi": "Valitse"
  },
  "stronglyLiked": {
    "eng": "Strongly liked",
    "fi": "Pidän paljon"
  },
  "veryUnfamiliar": {
    "eng": "Very unfamiliar",
    "fi": "Ei lainkaan tuttu"
  },
  "veryFamiliar": {
    "eng": "Very familiar",
    "fi": "Todella tuttu"
  },
  "thankYouRepeat": {
    "eng": "Well done! You reached midway! <br>Take a moment to move your wrists, take a breath, " +
    "and rest your eyes. <br>When you click 'Next', the task will continue.",
    "fi": "Kiitos vastauksista! Nyt tehtävä toistuu uudelleen."
  },
  "thankYouEnd": {
    "eng": "The task is now complete! Thank you for your participation",
    "fi": "Tehtävä on nyt valmis! Kiitos osallistumisesta."
  },
  "howOld": {
    "eng": "How old are you?",
    "fi": "Kuinka vanha olet?"
  },
  "gender": {
    "eng": "What is your gender?",
    "fi": "Mikä on sukupuolesi?"
  },
  "feedbackAsk": {
    "eng": "Do you have any feedback to give us about this task?",
    "fi": "Onko sinulla antaa meille palautetta tästä tehtävästä?"
  },
  "feedbackPlaceHolder": {
    "eng": "Was there any difficulty that you would like to report?",
    "fi": "Löysitkö ongelmia, joista haluaisit ilmoittaa meille?"
  },
  "genderOptions": [
    {
      "fi": "Nainen",
      "eng": "Woman"
    },
    {
      "fi": "Mies",
      "eng": "Man"
    },
    {
      "fi": "Joku muu",
      "eng": "Other"
    },
    {
      "fi": "En halua kertoa",
      "eng": "Prefer not to say"
    },
  ],
  "hearingLoss":
    {
      "fi": "Onko sinulla kuulonalenema?",
      "eng": "Do you have any hearing loss?"
    }
  ,
  "hearingOptions": [
    {
      "fi": "Ei kuulonalenemaa",
      "eng": "No hearing loss"
    },
    {
      "fi": "Lievä kuulonalenema",
      "eng": "Mild hearing loss"
    },
    {
      "fi": "Keskivaikea tai vaikea kuulonalenema",
      "eng": "Moderate to severe hearing loss"
    }
  ],
  "musicianship":{
    "eng":"Which title best describes you?",
    "fi":"Mikä seuraavista kuvaa sinua parhaiten?"
  },
  "musicianshipOptions":[
    {
      "fi":"Ei-muusikko",
      "eng":"Nonmusician"
    },
    {
      "fi":"Musiikkia rakastava ei-muusikko",
      "eng":"Music-loving nonmusician"
    },
    {
      "fi":"Harrastajamuusikko",
      "eng":"Amateur musician"
    },
    {
      "fi":"Aktiivinen harrastajamuusikko",
      "eng":"Serious amateur musician"
    },
    {
      "fi":"Puoliammattilaismuusikko",
      "eng":"Semiprofessional musician"
    },
    {
      "fi":"Ammattilaismuusikko",
      "eng":"Professional musician"
    }
  ],
  "education":{
    "fi":"Mikä on korkein kolulutustasosi?",
    "eng":"What is your highest level of education completed?"
  },
  "educationOptions":[
    {
      "fi":"Peruskoulu",
      "eng":"Middle school"
    },
    {
      "fi":"Lukio",
      "eng":"High school"
    },
    {
      "fi":"Ammattikoulu",
      "eng":"College / vocational training"
    },
    {
      "fi":"Alempi korkeakoulututkinto",
      "eng":"Bachelor's degree"
    },
    {
      "fi":"Ylempi korkeakoulututkinto",
      "eng":"Master's degree"
    },
    {
      "fi":"Tohtorintutkinto",
      "eng":"Doctoral degree"
    },
    {
      "fi":"Muu",
      "eng":"Other / prefer not to say"
    }
  ]

}
var consent_form = [
  "<p style='font-size:20px;text-align:left;line-height: normal'> MUSICAL EMOTIONS – PILOT STUDY </p>" +
    "<p style='font-size:14px;text-align:left;line-height: normal'>We ask you to participate in “Musical emotions – pilot study”, which investigates how suitable " + 
    "the items of a novel test are for measuring the perception of emotions in music. You are invited to the study " +
    "because you are an adult over 18 years old with normal hearing and with access to the internet, a device such as a " +
    "smartphone or laptop, and headphones. The study will involve around 150 adult participants."+
    "<br><br>This research notification describes the study and related participation. "+
    "You can find more information on the processing of your personal data in the following link: "+
    "https://jyu-my.sharepoint.com/:b:/g/personal/maarhart_jyu_fi/EcnAVidj975ImyASoMlEP4wBXerhPQ4FJCFAeXNpfYpUPw?e=f5h0uA"+
    "<br><t style='font-size:20px;text-align:left'> <br>1. Participation </t><br>Participation in this study is voluntary. You can refuse to participate in the study, stop participating "+
    "or cancel your previously given consent, without stating any reason for this and at any time during the study. This will have no negative "+
    "consequences to you. The study takes place during a lecture, but no registry will be made concerning which students were willing/refused to "+
    "participate. In other words, your willingness/refusal to take part of this study will not affect you or your situation in this course in any way. "+ 
    "If you stop participating in the study or if you cancel your consent, the personal data, samples and other information collected on you up to that point "+
    "will be used as part of the research material as far as it is necessary in order to ensure relevant research outcomes. "+  
    "<br><br><t style='font-size:20px;text-align:left'>2. Progress of the study </t><br>In this study, you will hear short clips of music and be asked which emotion it represents the most. To answer, you "+
    "will select one emoticon, with an emotion word written under it (eg. Angry, sad, happy, tender). Each round has 30 clips for you to rate and lasts around "+
    "10 minutes. The study has two rounds, so it takes around 20 minutes in total. Before the music listening part starts, you will also be asked some basic "+
    "demographics through a short survey (e.g. age and gender).<br><br><t style='font-size:20px;text-align:left'> 3. Possible benefits from the study</t><br>Taking part of the study might not bring any personal "+
    "benefit. The results of the study will benefit science, in the sense that it will help further develop an instrument to measure perception of emotions in music. "+
    "<br><br><t style='font-size:20px;text-align:left'>4. Possible risks and harm</t><br>Participation in the study is not expected to cause any risks, "+
    "harm or inconvenience. <br><br><t style='font-size:20px;text-align:left'>5. Study-related costs and compensations</t><br>No rewards will be paid for participation in the study. "+   
    "<br><br><t style='font-size:20px;text-align:left'>6. Informing about research results and research outcomes</t><br>The study will yield scientific publications, conference and seminar presentations, teaching material, and "+
    "practical applications (development of measurement instrument). No personal data will be published. Results will be treated in an aggregated manner, and subjects "+
    "cannot be identified from publications.<br><br><t style='font-size:20px;text-align:left'>7. Insurance coverage for research subjects</t><br>The University of Jyväskylä has insurances for its activities and research subjects. "+
    "The set of insurance includes a malpractice insurance, an operational liability insurance, and an optional insurance against accidents. During the research activities, the "+
    "subjects are covered by the insurance for accidents, damages and injuries inflicted by an external cause. The accident insurance is valid during measurements and on trips integrally "+
    "connected to them.<br><br><t style='font-size:20px;text-align:left'>8. Contact person for further information</t><br>Petri Toiviainen, petri.toiviainen@jyu.fi, +358503541753<br>Department of Music, Arts and Culture Studies, PO Box "+
    "35(M), 40014 University of Jyväskylä, Finland</p>",
    "<p style='font-size:20px;text-align:left;line-height: normal'> CONSENT TO PARTICIPATE IN SCIENTIFIC RESEARCH <br><br>" +
    "Musical emotions – pilot study </p><p style='font-size:14px;text-align:left;line-height: normal'> I understand that participation in the study is voluntary and that I can stop participating at any time, without giving a reason. " +
    "There will be no negative consequences for me if I withdraw. The data collected about me up to the point of withdrawal may still be used in the study. " +
    "<br><br>I have been adequately informed about the study and the processing of my personal data. I have received the information sheet about the study,  " +
    "as well as the privacy notice. I have also had the opportunity to ask the researchers further questions. " +
    "<br><br> By clicking the button, <br>- I accept that data will be collected from me as described in information sheet, " +
    "<br>   - I accept that my data is used in accordance with the procedures outlined in the privacy notice, " + 
    "<br>   - I confirm that I understand the information that I have received, " +
    "<br>   - I agree to participate in this study. </p> "
]
