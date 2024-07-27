//-------------------------------------------------------------------- JENKS NATURAL BREAKS ALGORITHM ---------------------------------------------------------
// # [Jenks natural breaks optimization](http://en.wikipedia.org/wiki/Jenks_natural_breaks_optimization)
//
// Implementations: [1](http://danieljlewis.org/files/2010/06/Jenks.pdf) (python),
// [2](https://github.com/vvoovv/djeo-jenks/blob/master/main.js) (buggy),
// [3](https://github.com/simogeo/geostats/blob/master/lib/geostats.js#L407) (works)
function jenks(data, n_classes) {
  // Compute the matrices required for Jenks breaks. These matrices
  // can be used for any classing of data with `classes <= n_classes`
  function getMatrices(data, n_classes) {
    // in the original implementation, these matrices are referred to
    // as `LC` and `OP`
    //
    // * lower_class_limits (LC): optimal lower class limits
    // * variance_combinations (OP): optimal variance combinations for all classes
    var lower_class_limits = [],
      variance_combinations = [],
      // loop counters
      i,
      j,
      // the variance, as computed at each step in the calculation
      variance = 0;

    // Initialize and fill each matrix with zeroes
    for (i = 0; i < data.length + 1; i++) {
      var tmp1 = [],
        tmp2 = [];
      for (j = 0; j < n_classes + 1; j++) {
        tmp1.push(0);
        tmp2.push(0);
      }
      lower_class_limits.push(tmp1);
      variance_combinations.push(tmp2);
    }

    for (i = 1; i < n_classes + 1; i++) {
      lower_class_limits[1][i] = 1;
      variance_combinations[1][i] = 0;
      // in the original implementation, 9999999 is used but
      // since Javascript has `Infinity`, we use that.
      for (j = 2; j < data.length + 1; j++) {
        variance_combinations[j][i] = Infinity;
      }
    }

    for (var l = 2; l < data.length + 1; l++) {
      // `SZ` originally. this is the sum of the values seen thus
      // far when calculating variance.
      var sum = 0,
        // `ZSQ` originally. the sum of squares of values seen
        // thus far
        sum_squares = 0,
        // `WT` originally. This is the number of
        w = 0,
        // `IV` originally
        i4 = 0;

      // in several instances, you could say `Math.pow(x, 2)`
      // instead of `x * x`, but this is slower in some browsers
      // introduces an unnecessary concept.
      for (var m = 1; m < l + 1; m++) {
        // `III` originally
        var lower_class_limit = l - m + 1,
          val = data[lower_class_limit - 1];

        // here we're estimating variance for each potential classing
        // of the data, for each potential number of classes. `w`
        // is the number of data points considered so far.
        w++;

        // increase the current sum and sum-of-squares
        sum += val;
        sum_squares += val * val;

        // the variance at this point in the sequence is the difference
        // between the sum of squares and the total x 2, over the number
        // of samples.
        variance = sum_squares - (sum * sum) / w;

        i4 = lower_class_limit - 1;

        if (i4 !== 0) {
          for (j = 2; j < n_classes + 1; j++) {
            // if adding this element to an existing class
            // will increase its variance beyond the limit, break
            // the class at this point, setting the lower_class_limit
            // at this point.
            if (
              variance_combinations[l][j] >=
              variance + variance_combinations[i4][j - 1]
            ) {
              lower_class_limits[l][j] = lower_class_limit;
              variance_combinations[l][j] =
                variance + variance_combinations[i4][j - 1];
            }
          }
        }
      }

      lower_class_limits[l][1] = 1;
      variance_combinations[l][1] = variance;
    }

    // return the two matrices. for just providing breaks, only
    // `lower_class_limits` is needed, but variances can be useful to
    // evaluage goodness of fit.
    return {
      lower_class_limits: lower_class_limits,
      variance_combinations: variance_combinations,
    };
  }

  // the second part of the jenks recipe: take the calculated matrices
  // and derive an array of n breaks.
  function breaks(data, lower_class_limits, n_classes) {
    var k = data.length - 1,
      kclass = [],
      countNum = n_classes;

    // the calculation of classes will never include the upper and
    // lower bounds, so we need to explicitly set them
    kclass[n_classes] = data[data.length - 1];
    kclass[0] = data[0];

    // the lower_class_limits matrix is used as indexes into itself
    // here: the `k` variable is reused in each iteration.
    while (countNum > 1) {
      kclass[countNum - 1] = data[lower_class_limits[k][countNum] - 2];
      k = lower_class_limits[k][countNum] - 1;
      countNum--;
    }

    return kclass;
  }

  if (n_classes > data.length) return null;

  // sort data in numerical order, since this is expected
  // by the matrices function
  data = data.slice().sort(function (a, b) {
    return a - b;
  });

  // get our basic matrices
  var matrices = getMatrices(data, n_classes),
    // we only need lower class limits here
    lower_class_limits = matrices.lower_class_limits;

  // extract n_classes out of the computed matrices
  return breaks(data, lower_class_limits, n_classes);
}

//-------------------------------------------------------------------- CONSTANTS ---------------------------------------------------------
const BREAKS_AMOUNT = 6;

const healthMetrics = [
  "Lack of health insurance crude prevalence (%)",
  "Binge drinking crude prevalence (%)",
  "Current smoking crude prevalence (%)",
  "Physical inactivity crude prevalence (%)",
  "Sleep <7 hours crude prevalence (%)",
  "Current asthma crude prevalence (%)",
  "High blood pressure crude prevalence (%)",
  "Cancer (except skin) crude prevalence (%)",
  "Cholesterol screening crude prevalence (%)",
  "Chronic kidney disease crude prevalence (%)",
  "Arthritis crude prevalence (%)",
  "Coronary heart disease crude prevalence (%)",
  "Diabetes crude prevalence (%)",
  "Obesity crude prevalence (%)",
  "Stroke crude prevalence (%)",
  "Annual checkup crude prevalence (%)",
  "Dental visit crude prevalence (%)",
  "Mammography use crude prevalence (%)",
  "Cervical cancer screening crude prevalence (%)",
  "Colorectal cancer screening crude prevalence (%)",
  "Depression crude prevalence (%)",
  "Frequent mental health distress crude prevalence (%)",
  "Frequent physical health distress crude prevalence (%)",
  "Fair or poor health crude prevalence (%)",
  "Any disability crude prevalence (%)",
  "Hearing disability crude prevalence (%)",
  "Vision disability crude prevalence (%)",
  "Cognitive disability crude prevalence (%)",
  "Mobility disability crude prevalence (%)",
  "Self-care disability crude prevalence (%)",
  "Independent living disability crude prevalence (%)",
];

const sunsetParkCensusTracts = [
  "36047000200",
  "36047001801",
  "36047001802",
  "36047001803",
  "36047001804",
  "36047002000",
  "36047002200",
  "36047007200",
  "36047007400",
  "36047007600",
  "36047007800",
  "36047008000",
  "36047008200",
  "36047008400",
  "36047010100",
  "36047014300",
  "36047014500",
  "36047014700",
  "36047008600",
  "36047008800",
  "36047009001",
  "36047009201",
  "36047009401",
  "36047009600",
  "36047009800",
  "36047010000",
  "36047010200",
  "36047010401",
  "36047010601",
  "36047010801",
  "36047011800",
  "36047012200",
  "36047009002",
  "36047009202",
  "36047009402",
  "36047010402",
  "36047010602",
  "36047010802",
  "36047011000",
  "36047011200",
  "36047011400",
  "36047011600",
];
const sunsetParkCensusTractsSet = new Set(sunsetParkCensusTracts);

const languageBreaks = {};
const healthData = {};
const healthBreaks = {};

const boroughs = ["Brooklyn", "Manhattan", "Queens", "Staten Island", "Bronx"];

const boroughLanguageData = {
  Brooklyn: {},
  Manhattan: {},
  Queens: {},
  "Staten Island": {},
  Bronx: {},
};

const boroughHealthData = {
  Brooklyn: {},
  Manhattan: {},
  Queens: {},
  "Staten Island": {},
  Bronx: {},
};

const boroughLanguageBreaks = {
  Brooklyn: {},
  Manhattan: {},
  Queens: {},
  "Staten Island": {},
  Bronx: {},
};

const boroughHealthBreaks = {
  Brooklyn: {},
  Manhattan: {},
  Queens: {},
  "Staten Island": {},
  Bronx: {},
};

const languageData = {
  Arabic: [],
  Chinese: [],
  French: [],
  German: [],
  Korean: [],
  Other: [],
  Other_Asia: [],
  Other_Indo: [],
  Russian: [],
  Spanish: [],
  Tagalog: [],
  Vietnamese: [],
  Total_pop: [],
};

const sunsetParkLanguageData = JSON.parse(JSON.stringify(languageData));
const sunsetParkHealthData = {};
const sunsetParkLanguageBreaks = {};
const sunsetParkHealthBreaks = {};

boroughs.forEach((borough) => {
  boroughLanguageData[borough] = JSON.parse(JSON.stringify(languageData));
  healthMetrics.forEach((metric) => {
    boroughHealthData[borough][metric] = [];
    sunsetParkHealthData[metric] = [];
  });
});

//---------------------------------------------------- Language and Demographics -----------------------------------------

languageGeoJsonData.features.forEach((feature) => {
  var properties = feature.properties;

  Object.keys(languageData).forEach((language) => {
    var value = parseFloat(properties[language]);
    if (!isNaN(value)) {
      languageData[language].push(value);

      boroughs.forEach((borough) => {
        if (properties.boroname.trim() === borough) {
          boroughLanguageData[borough][language].push(value);
        }
      });

      if (sunsetParkCensusTractsSet.has(String(properties.Tract))) {
        sunsetParkLanguageData[language].push(value);
      }
    }
  });
});

Object.keys(languageData).forEach((language) => {
  languageBreaks[language] = jenks(languageData[language], BREAKS_AMOUNT);
});

boroughs.forEach((borough) => {
  Object.keys(boroughLanguageData[borough]).forEach((language) => {
    boroughLanguageBreaks[borough][language] = jenks(
      boroughLanguageData[borough][language],
      BREAKS_AMOUNT
    );
  });
});

Object.keys(sunsetParkLanguageData).forEach((language) => {
  sunsetParkLanguageBreaks[language] = jenks(
    sunsetParkLanguageData[language],
    BREAKS_AMOUNT
  );
});

//---------------------------------------------------- Health -----------------------------------------

healthMetrics.forEach((metric) => (healthData[metric] = []));

healthDataGeojson.features.forEach((feature) => {
  var properties = feature.properties;

  healthMetrics.forEach((metric) => {
    var value = parseFloat(properties[metric]);
    if (!isNaN(value)) {
      healthData[metric].push(value);

      boroughs.forEach((borough) => {
        let countyNameMap = {
          Brooklyn: "Kings",
          Manhattan: "New York",
          Queens: "Queens",
          "Staten Island": "Richmond",
          Bronx: "Bronx",
        };

        if (properties["County name"] === countyNameMap[borough]) {
          boroughHealthData[borough][metric].push(value);
        }
      });

      if (sunsetParkCensusTractsSet.has(properties["Census tract FIPS"])) {
        sunsetParkHealthData[metric].push(value);
      }
    }
  });
});

Object.keys(healthData).forEach((metric) => {
  healthBreaks[metric] = jenks(healthData[metric], BREAKS_AMOUNT);
});

boroughs.forEach((borough) => {
  Object.keys(boroughHealthData[borough]).forEach((metric) => {
    boroughHealthBreaks[borough][metric] = jenks(
      boroughHealthData[borough][metric],
      BREAKS_AMOUNT
    );
  });
});

Object.keys(sunsetParkHealthData).forEach((metric) => {
  sunsetParkHealthBreaks[metric] = jenks(
    sunsetParkHealthData[metric],
    BREAKS_AMOUNT
  );
});

//---------------------------------------------------- Display -----------------------------------------
var languageBreaksDiv = document.getElementById("language-breaks");
var healthBreaksDiv = document.getElementById("health-breaks");
var brooklynLanguageBreaksDiv = document.getElementById(
  "brooklyn-language-breaks"
);
var brooklynHealthBreaksDiv = document.getElementById("brooklyn-health-breaks");
var manhattanLanguageBreaksDiv = document.getElementById(
  "manhattan-language-breaks"
);
var manhattanHealthBreaksDiv = document.getElementById(
  "manhattan-health-breaks"
);
var queensLanguageBreaksDiv = document.getElementById("queens-language-breaks");
var queensHealthBreaksDiv = document.getElementById("queens-health-breaks");
var statenLanguageBreaksDiv = document.getElementById("staten-language-breaks");
var statenHealthBreaksDiv = document.getElementById("staten-health-breaks");
var bronxLanguageBreaksDiv = document.getElementById("bronx-language-breaks");
var bronxHealthBreaksDiv = document.getElementById("bronx-health-breaks");
var sunsetParkLanguageBreaksDiv = document.getElementById(
  "sunset-park-language-breaks"
);
var sunsetParkHealthBreaksDiv = document.getElementById(
  "sunset-park-health-breaks"
);

const boroughDivs = {
  Brooklyn: {
    language: brooklynLanguageBreaksDiv,
    health: brooklynHealthBreaksDiv,
  },
  Manhattan: {
    language: manhattanLanguageBreaksDiv,
    health: manhattanHealthBreaksDiv,
  },
  Queens: {
    language: queensLanguageBreaksDiv,
    health: queensHealthBreaksDiv,
  },
  "Staten Island": {
    language: statenLanguageBreaksDiv,
    health: statenHealthBreaksDiv,
  },
  Bronx: {
    language: bronxLanguageBreaksDiv,
    health: bronxHealthBreaksDiv,
  },
  "Sunset Park": {
    language: sunsetParkLanguageBreaksDiv,
    health: sunsetParkHealthBreaksDiv,
  },
};

Object.keys(languageBreaks).forEach((language) => {
  var h3 = document.createElement("h3");
  h3.textContent = language;
  languageBreaksDiv.appendChild(h3);

  var p = document.createElement("p");
  p.textContent = languageBreaks[language].join(", ");
  languageBreaksDiv.appendChild(p);
});

boroughs.forEach((borough) => {
  Object.keys(boroughLanguageBreaks[borough]).forEach((language) => {
    var div = boroughDivs[borough].language;
    var h3 = document.createElement("h3");
    h3.textContent = language;
    div.appendChild(h3);

    var p = document.createElement("p");
    p.textContent = boroughLanguageBreaks[borough][language].join(", ");
    div.appendChild(p);
  });

  Object.keys(boroughHealthBreaks[borough]).forEach((metric) => {
    var div = boroughDivs[borough].health;
    var h3 = document.createElement("h3");
    h3.textContent = metric;
    div.appendChild(h3);

    var p = document.createElement("p");
    p.textContent = boroughHealthBreaks[borough][metric].join(", ");
    div.appendChild(p);
  });
});

Object.keys(healthBreaks).forEach((metric) => {
  var h3 = document.createElement("h3");
  h3.textContent = metric;
  healthBreaksDiv.appendChild(h3);

  var p = document.createElement("p");
  p.textContent = healthBreaks[metric].join(", ");
  healthBreaksDiv.appendChild(p);
});

Object.keys(sunsetParkLanguageBreaks).forEach((language) => {
  var h3 = document.createElement("h3");
  h3.textContent = language;
  sunsetParkLanguageBreaksDiv.appendChild(h3);

  var p = document.createElement("p");
  p.textContent = sunsetParkLanguageBreaks[language].join(", ");
  sunsetParkLanguageBreaksDiv.appendChild(p);
});

Object.keys(sunsetParkHealthBreaks).forEach((metric) => {
  var h3 = document.createElement("h3");
  h3.textContent = metric;
  sunsetParkHealthBreaksDiv.appendChild(h3);

  var p = document.createElement("p");
  p.textContent = sunsetParkHealthBreaks[metric].join(", ");
  sunsetParkHealthBreaksDiv.appendChild(p);
});
