  function line_pixel_length(d, theta, n) {
    // for angle in [pi/4,3*pi/4],
    // flip along diagonal (transpose) and call recursively
    if(theta > Math.PI /4 && theta < 3/4*Math.PI) {
      return line_pixel_length(d, Math.PI/2-theta, n);
    }

    // for angle in [3*pi/4,pi],
    // redefine line to go in opposite direction
    if(theta > Math.PI/2) {
      d =-d;
      theta = theta-Math.PI;
    }

    // for angle in [-Math.PI/4,0],
    //  flip along x-axis (up/down) and call recursively
    if(theta < 0) {
      // flip up down with reverse
      return flipud(line_pixel_length(-d,-theta,n));
    }

    if(theta > Math.PI/2 || theta < 0) {
      console.log("invalid angle")
      return;
    }

    var L = zeros(n,n);

    var ct = Math.cos(theta);
    var st = Math.sin(theta);

    var x0 = n/2 - d*st;
    var y0 = n/2 + d*ct;

    var y = y0 - x0*st/ct;
    var jy = Math.ceil(y);
    var dy = y+n % 1;

    for(var jx=1; jx <= n; jx++){
      var dynext = dy + st/ct;
      if(dynext < 1){
        if (jy >= 1 && jy <= n){
            L[n-jy][jx-1] = 1/ct;
        }
        dy = dynext;
      }else{
        if (jy>=1 && jy<=n) L[n-jy][jx-1] = (1 - dy)/st ;
        if (jy+1>=1 && jy+1<=n) L[n-(jy+1)][jx-1] = (dynext - 1)/st; 
        dy = dynext - 1;
        jy = jy + 1;
      }
    }
    return L;
  }

  function generateMeasurements(x, width, Nd, Ntheta, sigma) {
    var N = Nd*Ntheta;
    var y=zeros(N,1)           // will store the N measurements
    var lines_d=zeros(N,1);     // will store the position of each line
    var lines_theta=zeros(N,1); // will store the angle of each line
    var i=1;
    for(var itheta=1; itheta <= Ntheta; itheta++){
      for(var id=1; id <= Nd; id++){
        // equally spaced parallel lines, distance from first to
        // last is about 1.4*n_pixels (to ensure coverage of whole
        // image when theta=pi/4)
        lines_d[i-1]= 0.7*width*(id-Nd/2-0.5)/(Nd/2);
        
        // equally spaced angles from 0 to pi
        lines_theta[i-1]=Math.PI*itheta/Ntheta;

        // L is a matrix of the same size as the image
        // with entries giving the length of the line over each pixel
        var L = line_pixel_length(lines_d[i-1],lines_theta[i-1], width);
                 
        var l = col_flatten(L); // make matrix L into a vector, as for X
        y[i-1]= numeric.dot(l,x); //+ randn()*sigma
                 // l'*x gives "line integral" of line over image,
                 // that is, the intensity of each pixel is multiplied by the
                 // length of line over that pixel, and then add for all pixels;
                 // a random, Gaussian noise, with std sigma is added to the
                 // measurement
        i=i+1; // for next measurement line
      }
    }
    return {y: y, lines_d: lines_d, lines_theta: lines_theta}
  }

  function generateCoeffMatrix(x, n_pixels,  Nd, Ntheta, sigma) {
    var N = Nd*Ntheta;
    var meas = generateMeasurements(x,n_pixels,Nd,Ntheta,sigma);
    var A = zeros(N,Math.pow(n_pixels,2));
    for(var i = 1; i <= N; i++){
      var d = meas.lines_d[i-1];
      var theta = meas.lines_theta[i-1];
      var lengths = line_pixel_length(d, theta, n_pixels);
      lengths = col_flatten(lengths);
      A[i-1] = lengths;
    }
    return {A: numeric.ccsSparse(A), y: meas.y};
  }

  function genSmoothingMatrix(n_pixels) {
    var horzd = add(horzcat(eye(n_pixels-1), zeros(n_pixels-1,1)),
                    horzcat(zeros(n_pixels-1,1), timesScalar(eye(n_pixels-1),-1.0)));
    horzd = blkdiag(horzd,n_pixels);
    var g = n_pixels*(n_pixels-1);
    var vertd = add(horzcat(eye(g), zeros(g,n_pixels)),
                    horzcat(zeros(g,n_pixels), timesScalar(eye(g),-1.0)));
    var D = vertcat(horzd,vertd);
    return numeric.ccsSparse(D);
  }

  function matrixFromImage(name){
    // Get a reference to the image you want the pixels of and its dimensions
    var myImage = document.getElementById(name);
    var w = myImage.width, h = myImage.height;
    // Create a Canvas element
    var canvas = document.createElement('canvas');
    canvas.width = w;
    canvas.height = h;
    // Draw image onto the canvas
    var ctx = canvas.getContext('2d');
    ctx.drawImage(myImage, 0, 0);
    // Finally, get the image data
    // ('data' is an array of RGBA pixel values for each pixel)
    var idata = ctx.getImageData(0, 0, w, h);
    return matrix(h,w, function(i,j) {
      var idx = (i*w + j)*4;
      return (idata.data[idx]*0.2989 + 
             idata.data[idx+1]*0.5870 + 
             idata.data[idx+2]*0.1140); 
    });
  }

  function imageFromMatrix(X){
    var canvas = document.createElement("canvas");
    var ctx = canvas.getContext("2d");
    var w = X[0].length; var h = X.length;
    canvas.width = w;
    canvas.height = h;
    var imgData = ctx.getImageData(0,0,w,h);
    for(var y = 0; y < h; y++) {
      for(var x = 0; x < w; x++) {
        var pos = (y*w+ x) * 4; 
        imgData.data[pos  ] = X[y][x]; // R Value
        imgData.data[pos+1] = X[y][x]; // G value
        imgData.data[pos+2] = X[y][x]; // B value
        imgData.data[pos+3] = 255;     // alpha channel
      }
    }
    ctx.putImageData(imgData,0,0);
    var image = new Image();
    image.src = canvas.toDataURL();
    return image;
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
    var D = genSmoothingMatrix(n_pixels)
    D = numeric.dot(numeric.transpose(D),D);
    console.log("Generated smoothing matrix");
    var lambdas = numeric.pow(10,numeric.linspace(-5,5,10));
    for(i=0; i< lambdas.length; i++){
      var lambda = lambdas[i];
      var M = numeric.add(numeric.dot(numeric.transpose(sys.A),sys.A), numeric.mul(D,lambda));
      console.log("created final matrix");
      var b = numeric.dot(numeric.transpose(sys.A),sys.y);
      console.log("solving");
      var xhat = numeric.solve(M,b);
      console.log("solved");
      var Xhat = col_reshape(xhat);
      var img = imageFromMatrix(Xhat);
      console.log("adding image");
      document.getElementById("tomo-vis").appendChild(img);
    }
  }

