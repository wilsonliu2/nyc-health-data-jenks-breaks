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

//-------------------------------------------------------------------- DATAS ---------------------------------------------------------

/*
Arabic
Chinese
French
German
Korean
Other
Other_Asia
Other_Indo
Russian
Spanish
Tagalog
Vietnamese
Total_pop
Lack of health insurance crude prevalence (%)
Binge drinking crude prevalence (%)
Current smoking crude prevalence (%)
Physical inactivity crude prevalence (%)
Sleep <7 hours crude prevalence (%)
Current asthma crude prevalence (%)
High blood pressure crude prevalence (%)
Cancer (except skin) crude prevalence (%)
Cholesterol screening crude prevalence (%)
Chronic kidney disease crude prevalence (%)
Arthritis crude prevalence (%)
Coronary heart disease crude prevalence (%)
Diabetes crude prevalence (%)
Obesity crude prevalence (%)
Stroke crude prevalence (%)
Current asthma crude prevalence (%)
High blood pressure crude prevalence (%)
Cancer (except skin) crude prevalence (%)
Cholesterol screening crude prevalence (%)
Chronic kidney disease crude prevalence (%)
Arthritis crude prevalence (%)
Coronary heart disease crude prevalence (%)
Diabetes crude prevalence (%)
Obesity crude prevalence (%)
Stroke crude prevalence (%)
Annual checkup crude prevalence (%)
Dental visit crude prevalence (%)
Cholesterol screening crude prevalence (%)
Mammography use crude prevalence (%)
Cervical cancer screening crude prevalence (%)
Colorectal cancer screening crude prevalence (%)
Depression crude prevalence (%)
Frequent mental health distress crude prevalence (%)
Frequent physical health distress crude prevalence (%)
Fair or poor health crude prevalence (%)
Any disability crude prevalence (%)
Any disability crude prevalence (%)
Hearing disability crude prevalence (%)
Vision disability crude prevalence (%)
Cognitive disability crude prevalence (%)
Mobility disability crude prevalence (%)
Self-care disability crude prevalence (%)
Independent living disability crude prevalence (%)
*/

const BREAKS_AMOUNT = 6;

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

//---------------------------------------------------- Language and Demographics -----------------------------------------

languageGeoJsonData.features.forEach((feature) => {
  var properties = feature.properties;

  Object.keys(languageData).forEach((language) => {
    var value = parseFloat(properties[language]);
    if (!isNaN(value)) {
      languageData[language].push(value);
    }
  });
});

// All language and demographics data
// console.log(languageData);

var languageBreaks = {};

Object.keys(languageData).forEach((language) => {
  languageBreaks[language] = jenks(languageData[language], BREAKS_AMOUNT);
});

// Display language breaks
// console.log(languageBreaks);

//---------------------------------------------------- Health -----------------------------------------

var healthData = {};
healthMetrics.forEach((metric) => (healthData[metric] = []));

healthDataGeojson.features.forEach((feature) => {
  var properties = feature.properties;

  healthMetrics.forEach((metric) => {
    var value = parseFloat(properties[metric]);
    if (!isNaN(value)) {
      healthData[metric].push(value);
    }
  });
});

// All health data
// console.log(healthData);

var healthBreaks = {};

// Calculate Jenks breaks for each health metric
Object.keys(healthData).forEach((metric) => {
  healthBreaks[metric] = jenks(healthData[metric], BREAKS_AMOUNT);
});

// Display health breaks
// console.log(healthBreaks);

//---------------------------------------------------- Display -----------------------------------------
var languageBreaksDiv = document.getElementById("language-breaks");
var healthBreaksDiv = document.getElementById("health-breaks");

Object.keys(languageBreaks).forEach((language) => {
  var h3 = document.createElement("h3");
  h3.textContent = language;
  languageBreaksDiv.appendChild(h3);

  var p = document.createElement("p");
  p.textContent = languageBreaks[language].join(", ");
  languageBreaksDiv.appendChild(p);
});

Object.keys(healthBreaks).forEach((metric) => {
  var h3 = document.createElement("h3");
  h3.textContent = metric;
  healthBreaksDiv.appendChild(h3);

  var p = document.createElement("p");
  p.textContent = healthBreaks[metric].join(", ");
  healthBreaksDiv.appendChild(p);
});
