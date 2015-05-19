$(function() {
  var $rangeSlider = $('#range-slider');
 
  var numLineIntegrals = 25;
 
  rangeSlider($rangeSlider[0], {
      value: numLineIntegrals,
      drag: function(value) {
        // Update the randomness after the user drags the slider
        // and reset the points to be clustered
        numLineIntegrals= value;
        resetLineIntegrals();
      }
  });

  function resetLineIntegrals() {
    console.log(numLineIntegrals);

  }

  function tomo() {
    var Nd = 40;       // number of parallel lines for each angle
    var Ntheta = 40    //  number of angles (equally spaced from 0 to pi)
    var N = Nd*Ntheta  // total number of lines (i.e. number of measurements
    var sigma = 20     // noise level (standard deviation for normal dist.)
    var X = matrixFromImage("theimage");
    console.log("Loaded Image");
    var n_pixels = X.length;
    var x = col_flatten(X);
    var sys = generateCoeffMatrix(x, n_pixels, Nd, Ntheta, sigma);
    console.log("Generated Problem");
    var D = numeric.ccsSparse(genSmoothingMatrix(n_pixels))
    D = numeric.dot(numeric.transpose(D),D);
    console.log("Generated smoothing matrix");
    var lambdas = numeric.pow(10,numeric.linspace(-5,5,10));
    for(i=0; i< lambdas.length; i++){
      var lambda = lambdas[i];
      var M = numeric.add(numeric.dot(numeric.transpose(sys.A),sys.A), numeric.mul(D,lambda));
      console.log("created final matrix");
      var b = numeric.dot(numeric.transpose(sys.A),sys.y);
      console.log("solving");
      var xhat = solve(M,b);
      console.log("solved");
      var Xhat = col_reshape(xhat);
      var img = imageFromMatrix(Xhat);
      console.log("adding image");
      document.getElementById("tomo-vis").appendChild(img);
    }
  }

})
