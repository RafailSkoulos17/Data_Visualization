//Functions we implemented


// Function that computes the 1D cubic interpolation. g0,g1,g2,g3 contain the values of the voxels that we want to interpolate
    // factor contains the distance from the value g1 to the position we want to interpolate to.
    // We assume the out of bounce checks have been done earlier
    
    private float cubicinterpolate(float g0, float g1, float g2, float g3, float factor) {
            
        // compute the result of the 1D cubic interpolation
        float result = (float) ( a * ((-factor+2.f* Math.pow(factor, 2)-Math.pow(factor, 3))*g0
               +(2.f-5.f*Math.pow(factor, 2)+3.f*Math.pow(factor, 3))*g1
               +(factor+4.f* Math.pow(factor, 2)-3.f* Math.pow(factor, 3))*g2
               +(-Math.pow(factor, 2)+ Math.pow(factor, 3))*g3));

        return Math.abs(result); 
    }


    //--------------------------------------------------------------------------------------------


    // 2D cubic interpolation implemented here. We do it for plane XY. Coord contains the position.
    // We assume the out of bounce checks have been done earlier
    private float bicubicinterpolateXY(double[] coord, int z) {
            
        
        int x = (int) Math.floor(coord[0]);  
        int y = (int) Math.floor(coord[1]);
        
        //computethe factors for each axis
        float fac_x = (float) coord[0] - x;  
        float fac_y = (float) coord[1] - y;
        float fac_z = (float) coord[2] - z;

        // do the four cubic interpolations
        float t0 = cubicinterpolate(getVoxel(x-1, y-1, z), getVoxel(x, y-1, z), getVoxel(x+1, y-1, z), getVoxel(x+2, y-1, z), fac_x);
        float t1 = cubicinterpolate(getVoxel(x-1, y, z), getVoxel(x, y, z), getVoxel(x+1, y, z), getVoxel(x+2, y, z), fac_x);
        float t2 = cubicinterpolate(getVoxel(x-1, y+1, z), getVoxel(x, y+1, z), getVoxel(x+1, y+1, z), getVoxel(x+2, y+1, z), fac_x);
        float t3 = cubicinterpolate(getVoxel(x-1, y+2, z), getVoxel(x, y+2, z), getVoxel(x+1, y+2, z), getVoxel(x+2, y+2, z), fac_x);
        
        // compute the final result of the bicubic interpolation
        float result = cubicinterpolate(t0, t1, t2, t3, fac_y);
         
                 
        return result; 

    }

    //-----------------------------------------------------------------------------------


    // 3D cubic interpolation implemented here given a position in the volume given by coord.
    
    public float getVoxelTriCubicInterpolate(double[] coord) {
        if (coord[0] < 1 || coord[0] > (dimX-3) || coord[1] < 1 || coord[1] > (dimY-3)
                || coord[2] < 1 || coord[2] > (dimZ-3)) {
            return 0;
        }

        int z = (int) Math.floor(coord[2]);
        float fac_z = (float) coord[2] - z;
        
        // do the bicubic interpolations for z-1, z, z+1 and z+2
        float bicubic_0 = bicubicinterpolateXY(coord, z-1);
        float bicubic_1 = bicubicinterpolateXY(coord, z);
        float bicubic_2 = bicubicinterpolateXY(coord, z+1);
        float bicubic_3 = bicubicinterpolateXY(coord, z+2);
        
        // apply cubic interpolation on the results of the 4 previous bicubic interpolations
        // to get the final result of the tricubic interpolation
        float result = cubicinterpolate(bicubic_0, bicubic_1, bicubic_2, bicubic_3, fac_z);
              
        return result; 
        

    }

    // -----------------------------------------------------------------------------------------------------


    /This function linearly interpolates gradient vector g0 and g1 given the factor (t) 
//the result is given at result. You can use it to tri-linearly interpolate the gradient 
    
	private void interpolate(VoxelGradient g0, VoxelGradient g1, float factor, VoxelGradient result) {
            
        result.x = (1.0f - factor)*g0.x + factor*g1.x;
        result.y = (1.0f - factor)*g0.y + factor*g1.y;
        result.z = (1.0f - factor)*g0.z + factor*g1.z;

        result.mag = (float) Math.sqrt(result.x * result.x + result.y * result.y + result.z * result.z);
    }


    //--------------------------------------------------------------------------

    // This function should return linearly interpolated gradient for position coord[]
// right now it returns the nearest neighbour        
        
    public VoxelGradient getGradient(double[] coord) {
        
        if (coord[0] < 0 || coord[0] > (dimX-2) || coord[1] < 0 || coord[1] > (dimY-2)
                || coord[2] < 0 || coord[2] > (dimZ-2)) {
            return zero;
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);
        
        //find distance of interpolation vertex along the three axes
        float distx = (float) coord[0] - x;
        float disty = (float) coord[1] - y;
        float distz = (float) coord[2] - z;
        
        //get the gradient for all the vertices of interest
        VoxelGradient g000 = getGradient(x,y,z);
        VoxelGradient g001 = getGradient(x,y,z+1);
        VoxelGradient g010 = getGradient(x,y+1,z);
        VoxelGradient g011 = getGradient(x,y+1,z+1);
        VoxelGradient g100 = getGradient(x+1,y,z);
        VoxelGradient g101 = getGradient(x+1,y,z+1);
        VoxelGradient g110 = getGradient(x+1,y+1,z);
        VoxelGradient g111 = getGradient(x+1,y+1,z+1);
        
        //find interval gradient along axis x with linear interpolation
        VoxelGradient v0 = new VoxelGradient();
        VoxelGradient v1 = new VoxelGradient();
        VoxelGradient v2 = new VoxelGradient();
        VoxelGradient v3 = new VoxelGradient();
        
        //we do trilinear interpolation to calculate in the specific point (we just have the gradients of the vertices)
        interpolate(g000, g100, distx, v0);
        interpolate(g001, g101, distx, v1);
        interpolate(g010, g110, distx, v2);
        interpolate(g011, g111, distx, v3);
        
        //find interval vertices along axis y with linear interpolation
        VoxelGradient v4 = new VoxelGradient();
        VoxelGradient v5 = new VoxelGradient();
        interpolate(v0, v2, disty, v4);
        interpolate(v1, v3, disty, v5);
        
        
        //find interpolation vertex along z axis
        VoxelGradient v6 = new VoxelGradient();
        interpolate(v4, v5, distz, v6);
        //return the gradient of the requested sample
        return v6;

    }

    //---------------------------------------------------


    //Function that updates the "image" attribute using the Isosurface raycasting
    //It returns the color assigned to a ray/pixel given it's starting point (entryPoint) and the direction of the ray(rayVector).
    // exitPoint is the last point.
    //ray must be sampled with a distance defined by the sampleStep
   
   int traceRayIso(double[] entryPoint, double[] exitPoint, double[] rayVector, double sampleStep) {
       
        TFColor shading = new TFColor();
        TFColor newcolor = new TFColor();
        
        VoxelGradient gradient = new VoxelGradient();
        
        double[] viewVector = new double[3];
        VectorMath.setVector(viewVector, -rayVector[0], -rayVector[1], -rayVector[2]);
        double[] lightVector = new double[3];
        //We define the light vector as directed toward the view point (which is the source of the light)
        // another light vector would be possible
         VectorMath.setVector(lightVector, rayVector[0], rayVector[1], rayVector[2]);
         
        // defines if bisection will be used
        boolean bisection = false;
        //defines if there was a hit
        boolean found=false;
        
        //the current position is initialized as the entry point
        double[] currentPos = new double[3];
        VectorMath.setVector(currentPos, entryPoint[0], entryPoint[1], entryPoint[2]);
        double[] increments = new double[3];
        VectorMath.setVector(increments, rayVector[0] * sampleStep, rayVector[1] * sampleStep, rayVector[2] * sampleStep);
        double[] prevPos = new double[3];
        VectorMath.setVector(prevPos, currentPos[0], currentPos[1], currentPos[2]);
        double[] best_position = new double[3];
        
        float previousvalue = 0;
        float value = 0;
        
       // compute the number of times we need to sample
        int nrSamples = 1 + (int) Math.floor(VectorMath.distance(entryPoint, exitPoint) / sampleStep);
        
        //sample through the ray
        while (!found && nrSamples>0){
            previousvalue = value;
            //get the current intensity value
            value = volume.getVoxelLinearInterpolate(currentPos); 
            if (value >= iso_value) { //if there was a hit
                found=true; 
                if (bisection){ //if we use bisection mode
                    best_position = bisection_accuracy (currentPos, prevPos, previousvalue, value, iso_value);
                    gradient = gradients.getGradient(best_position);
                }
                else{
                    gradient = gradients.getGradient(currentPos); 
                }
                break;
            }
            // set the previous position
            VectorMath.setVector(prevPos, currentPos[0], currentPos[1], currentPos[2]);
            
            // move to the next position
            for (int i = 0; i < 3; i++) {
                currentPos[i] += increments[i];
            }
            nrSamples--;
        }

        //Initialization of the colors as floating point values
        double alpha;
        
        // if there was a miss, make the voxel transarent
        if (found) { 
            alpha = 1.0;
        } else {
            alpha = 0.0;
        }             
        
        //isoColor contains the isosurface color from the interface
        newcolor.r = isoColor.r;
        newcolor.g = isoColor.g;
        newcolor.b =  isoColor.b;
        
        //if shading mode is enabled
        if(shadingMode){
                //if the voxel is not transparent we apply shading using Phong's model
                if(alpha==1.0){
                    //calculate the color and the opacity after the application of Phong Shading
                    shading=computePhongShading(newcolor, gradient, lightVector, viewVector);
                    
                    newcolor.r=shading.r;
                    newcolor.g=shading.g;
                    newcolor.b=shading.b;
                    newcolor.a=shading.a;
                    alpha = shading.a;

                }
            }
        
        //compute the color
        int color = computeImageColor(newcolor.r,newcolor.g,newcolor.b,alpha);
        return color;
    }


    //----------------------------------------------------------------------------------


    // Given the current sample position, increment vector of the sample (vector from previous sample to current sample) and sample Step. 
   // Previous sample value and current sample value, isovalue value
   // The function should search for a position where the iso_value passes that it is more precise.
   // Some unused argumenets were removed
   double[] bisection_accuracy (double[] currentPos, double[] prevPos, float previousvalue,float value, float iso_value)
   {    
       
        double[] middlePos = new double[3];
        // define the threshold we want for the difference between iso_value and the current value
        double isoTol = 1;
        // define maximum number of bisections
        
        int numBisections = 20;
        int middleValue = 0;
       
       for (int i=0; i<numBisections; i++){
            
           // get the middle point between current and previous position
            VectorMath.setVector(middlePos, (prevPos[0] + currentPos[0])/2.0, (prevPos[1] + currentPos[1])/2.0, (prevPos[2] + currentPos[2])/2.0);
            // get the intensity for the middle point
            middleValue = (int) volume.getVoxelLinearInterpolate(middlePos);
            
            // compares the isovalue and the value in the middle point
            // and it devided the interval in the right way
            if (Math.abs(middleValue - iso_value) <= isoTol){ //if we are close to the isovalue return this position
                return middlePos;
            }
            else if (middleValue > iso_value){ // if the value in the middle point is bigger than isovalue, 
                if (previousvalue > value){ // if the intensity is decreasing
                    prevPos = middlePos;
                }
                else{
                    currentPos = middlePos;
                }
            }
            else if (middleValue < iso_value){ // if the value in the middle point is smaller than isovalue
                if (previousvalue > value){ // if the intensity is decreasing
                    currentPos = middlePos; 
                }
                else{
                     prevPos = middlePos;     
                }
            }     
        }
        return middlePos;
   }

   //--------------------------------------------------------------------------------------------


   //Function that updates the "image" attribute using the compositing// accumulatted raycasting
    //It returns the color assigned to a ray/pixel given it's starting point (entryPoint) and the direction of the ray(rayVector).
    // exitPoint is the last point.
    //ray must be sampled with a distance defined by the sampleStep
   
    TFColor BackToFront(TFColor AccumulatedColor, TFColor newColor, double alpha){
        AccumulatedColor.r= newColor.r * alpha + (1.0 - alpha) * AccumulatedColor.r;
        AccumulatedColor.g= newColor.g * alpha + (1.0 - alpha) * AccumulatedColor.g;
        AccumulatedColor.b= newColor.b * alpha + (1.0 - alpha) * AccumulatedColor.b;
        AccumulatedColor.a = alpha + (1.0 - alpha) * AccumulatedColor.a;
     return AccumulatedColor;
   }


   //-------------------------------------------------------------------------------------
   
    int traceRayComposite(double[] entryPoint, double[] exitPoint, double[] rayVector, double sampleStep) {
        
        TFColor voxel_color = new TFColor();
        TFColor newColor=new TFColor();
        TFColor shading=new TFColor();
        
        //initialize variable color where we store the result of color composition
        TFColor AccumulatedColor=new TFColor();
        
        //initialize the light vector
        double[] lightVector = new double[3];
        VectorMath.setVector(lightVector, rayVector[0], rayVector[1], rayVector[2]);
        
        //initialize the view vector
        double[] viewVector = new double[3];
        VectorMath.setVector(viewVector, -rayVector[0], -rayVector[1], -rayVector[2]);

        //compute the increment and the number of samples
        double[] increments = new double[3];
        VectorMath.setVector(increments, rayVector[0] * sampleStep, rayVector[1] * sampleStep, rayVector[2] * sampleStep);
        
        
        //the current position is initialized as the entry point
        double[] currentPos = new double[3];
        double voxel_intensity;
        double alpha = 0; 
        
        // Define the distance and direction between each sample step along the three axes
    	double stepDistX = -rayVector [0] * sampleStep;
        double stepDistY = -rayVector [1] * sampleStep;
        double stepDistZ = -rayVector [2] * sampleStep;

        // Define array that holds the 3 coordinates of each sample.
        double[] voxelCoord = new double[3];
        //define the first sample which will be the last one (due to back to front color composing)
        voxelCoord[0] = exitPoint[0];
        voxelCoord[1] = exitPoint[1];
        voxelCoord[2] = exitPoint[2];
        
        double[] prevCoord = new double[3];
        double[] nexCoord = new double[3];
        
        VectorMath.setVector(currentPos, exitPoint[0], exitPoint[1], exitPoint[2]);
        int resultColor=0; 

        // Compute the number of times we need to sample
        int nrSamples = 1 + (int) Math.floor(VectorMath.distance(entryPoint, exitPoint) / sampleStep);
        
        //initialize all colors and opacity with zero
        AccumulatedColor.r=0.0;
        AccumulatedColor.b=0.0;
        AccumulatedColor.g=0.0;
        AccumulatedColor.a=0.0;
        
        

        
        //the light vector is directed toward the view point (which is the source of the light)
        // itarate through the ray
        do {

            // get the gradient and the intensity of the voxel
            voxel_intensity = volume.getVoxelLinearInterpolate(currentPos);
            VoxelGradient gradient = gradients.getGradient(currentPos);

            if (compositingMode){
                //based on the sample intensity, we get its color and opacity through the transfer function
                voxel_color = tFunc.getColor((int)voxel_intensity);
                alpha = voxel_color.a;
                //assign the color to our variable newColor
                newColor.r =voxel_color.r;
                newColor.g =voxel_color.g;
                newColor.b =voxel_color.b;
                newColor.a =voxel_color.a;
            }
            
            //if 2d transfer function mode is enabled
            if (tf2dMode){
                
                //get the color of the voxel via the 2d transfer function
                voxel_color=tFunc2D.color;
                
                 //assign the color to our variable newColor
                newColor.a = voxel_color.a;
                newColor.r = voxel_color.r;
                newColor.g = voxel_color.g;
                newColor.b = voxel_color.b;
                
                //for this mode we have to calculate the opacity from the triangular classification widget
                //the following function gets the value of the material, the factor r, the voxel value and the gradient magnitude
                // and returns the opacity to be assgned
                alpha = computeOpacity2DTF(tFunc2D.baseIntensity,tFunc2D.radius, voxel_intensity, gradient.mag);
                newColor.a = alpha;
            }
            
            //if shading mode is enabled
            if(shadingMode){
                
                //if the voxel is not transparent we apply shading using Phong's shading model
                if(alpha>0.0){
                    //calculate the color and opacity after the application of phong shading
                    shading=computePhongShading(newColor, gradient, lightVector, viewVector);
                    
                    //assign the color to our variable newColor
                    newColor.r=shading.r;
                    newColor.g=shading.g;
                    newColor.b=shading.b;
                    newColor.a=shading.a;
                    alpha = shading.a;

                }
            }
            
            // apply back to front compositing
            AccumulatedColor=BackToFront(AccumulatedColor,newColor,alpha);
           
            //implementation of early ray termination to speed up our calculations
            if (newColor.a>=0.98 && newColor.a<=1.0){
            	//if opacity of sample is approximately 1.0 we stop calculating samples along the ray
                break;
            }
            //update current position
            for (int i = 0; i < 3; i++) {
                currentPos[i] -= increments[i]; 
            }
                        
            nrSamples--;
        } while (nrSamples > 0);
        
        AccumulatedColor.a = 1.0;
       //compute the final color to be assigned
       resultColor=computeImageColor(AccumulatedColor.r,AccumulatedColor.g,AccumulatedColor.b,AccumulatedColor.a);
       return resultColor;
    }

//---------------------------------------------------------------------------------

TFColor computePhongShading(TFColor voxel_color, VoxelGradient gradient, double[] lightVector,
         double[] viewVector) {
        
        // initialize the shading color for the given voxel
	TFColor shadingColor=new TFColor(voxel_color.r,voxel_color.g,voxel_color.b,voxel_color.a);
        
        //parameters as given by the exercise
        double kd=0.7;
	double ka=0.1;
	double ks=0.2;
	double n_power=100;

        //normal vector 
	double[] normal= new double[3];
        //normalize the gradient by dividing it by its magnitude and keep it to the normal vector
	VectorMath.setVector(normal,gradient.x/gradient.mag, gradient.y/gradient.mag, gradient.z/gradient.mag);

	//calculate ambient and diffuse
        
        //calculate the cosine of the angle formed by lightVector and normal vector
	double diffuse=VectorMath.dotproduct(normal,lightVector);

        // calculate the color that results from ambient and diffuse based on the equation Ld*kd*ci,d(L.n)+La*ka*ci,a
        // where we assumed that La,Ld are equal to one
	if (diffuse > 0){
            shadingColor.r = shadingColor.r * diffuse * kd + voxel_color.r * ka ;
            shadingColor.g = shadingColor.g * diffuse * kd + voxel_color.g * ka;
            shadingColor.b = shadingColor.b * diffuse * kd + voxel_color.b * ka;
	}

	//calculate the r vector based on the equation  R=2(n.L)*n-L
      	double [] rVector=new double[3];
        for(int i=0;i<3;i++){
            rVector[i]=2*diffuse*normal[i]-lightVector[i];
        }

	//Calculate specular and VR product(this product is equal to the cosine of the angle between the viewVector and r vector)
	double specular=VectorMath.dotproduct(normal,rVector);
        double vrProduct=VectorMath.dotproduct(viewVector,rVector);

        //Calculate specular based on the equation Ls*ks*ci,s*(V.R)^n
	if (specular > 0) {
            //following is the procedure for specular based on Phong's model
//            shadingColor.r += ks * Math.pow(vrProduct, n_power);
//            shadingColor.g += ks * Math.pow(vrProduct, n_power);
//            shadingColor.b += ks * Math.pow(vrProduct, n_power);
            
            //we tried implemented an alternative for Phong's Model which is faster
            shadingColor.r += ks * ((vrProduct/(n_power-n_power*vrProduct +vrProduct)));
            shadingColor.g += ks * ((vrProduct/(n_power-n_power*vrProduct +vrProduct)));
            shadingColor.b += ks * ((vrProduct/(n_power-n_power*vrProduct +vrProduct)));
            
            // we used the absolute values, because of some negative values which were observed
            shadingColor.r = Math.abs(shadingColor.r);
            shadingColor.g = Math.abs(shadingColor.g);
            shadingColor.b = Math.abs(shadingColor.b);

        }

        shadingColor.r = shadingColor.r > 1.0 ? 1.0 : shadingColor.r;
        shadingColor.g = shadingColor.g > 1.0 ? 1.0 : shadingColor.g;
        shadingColor.b = shadingColor.b > 1.0 ? 1.0 : shadingColor.b;
        
        return shadingColor;
    }
    

    // Implements the basic tracing of rays trough the image and given the
    // camera transformation
    // It calls the functions depending on the raycasting mode
  
    void raycast(double[] viewMatrix) {
        
    	//data allocation
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        double[] pixelCoord = new double[3];
        double[] entryPoint = new double[3];
        double[] exitPoint = new double[3];
        
        // increment in the pixel domain in pixel units
        int increment = 1;
        // sample step in voxel units
        int sampleStep = 1;
        
         //We increases the sample step to 5 so we get less samples per ray 
        //and the number of pixels we look at a time to 2, when user is interacting with the simulator
        //This way the simulator is more responsive and does not lag
        if (this.interactiveMode){
            increment = 2;
            sampleStep = 5;
        }
        
        // reset the image to black
        resetImage();
        

        // vector uVec and vVec define the view plane, 
        // perpendicular to the view vector viewVec which is going from the view point towards the object
        // uVec contains the up vector of the camera in world coordinates (image vertical)
        // vVec contains the horizontal vector in world coordinates (image horizontal)
        getViewPlaneVectors(viewMatrix,viewVec,uVec,vVec);
        
        
        // The result of the visualization is saved in an image(texture)
        // we update the vector according to the resolution factor
        // If the resolution is 0.25 we will sample 4 times more points. 
        for(int k=0;k<3;k++)
        {
            uVec[k]=res_factor*uVec[k];
            vVec[k]=res_factor*vVec[k];
        }
        
        
        
       // We get the size of the image/texture we will be puting the result of the 
        // volume rendering operation.
        int imageW=image.getWidth();
        int imageH=image.getHeight();

        int[] imageCenter = new int[2];
        // Center of the image/texture 
        imageCenter[0]= imageW/2;
        imageCenter[1]= imageH/2;
        
        // imageW/ image H contains the real width of the image we will use given the resolution. 
        //The resolution is generated once based on the maximum resolution.
        imageW = (int) (imageW*((max_res_factor/res_factor)));
        imageH = (int) (imageH*((max_res_factor/res_factor)));
    
        //The rayVector is pointing towards the scene
        double[] rayVector = new double[3];
        rayVector[0]=-viewVec[0];
        rayVector[1]=-viewVec[1];
        rayVector[2]=-viewVec[2];

             
        // compute the volume center
        double[] volumeCenter = new double[3];
        computeVolumeCenter(volumeCenter);

        
        // ray computation for each pixel
        for (int j = imageCenter[1] - imageH/2; j < imageCenter[1] + imageH/2; j += increment) {
            for (int i =  imageCenter[0] - imageW/2; i <imageCenter[0] + imageW/2; i += increment) {
                // compute starting points of rays in a plane shifted backwards to a position behind the data set
            	computePixelCoordinatesBehindFloat(pixelCoord,viewVec,uVec,vVec,i,j);
            	// compute the entry and exit point of the ray
                computeEntryAndExit(pixelCoord, rayVector, entryPoint, exitPoint);
                if ((entryPoint[0] > -1.0) && (exitPoint[0] > -1.0)) {
                    int val = 0;
                    if (compositingMode || tf2dMode) {
                        val = traceRayComposite(entryPoint, exitPoint, rayVector, sampleStep);
                    } else if (mipMode) {
                        val = traceRayMIP(entryPoint, exitPoint, rayVector, sampleStep);
                    } else if (isoMode){
                        val= traceRayIso(entryPoint,exitPoint,rayVector, sampleStep);
                    }
                    for (int ii = i; ii < i + increment; ii++) {
                        for (int jj = j; jj < j + increment; jj++) {
                            image.setRGB(ii, jj, val);
                        }
                    }
                }

            }
        }
    }


//-------------------------------------------------



// Compute the opacity based on the value of the pixel and the values of the
// triangle widget tFunc2D contains the values of the baseintensity and radius
// tFunc2D.baseIntensity, tFunc2D.radius they are in image intensity units

 
    
    
public double computeOpacity2DTF(double material_value, double material_r,
        double voxelValue, double gradMagnitude) {
    
    //initialize opacity to zero
    double opacity = 0.0;
    
    //get the coordinbates of the triangle's vertexes a, b and c
    double[] a = {material_value - material_r, 173.1};
    double[] b = {material_value + material_r, 173.1};
    double[] c = {material_value, 0};
    
    //get the coordinates of the center of mass(centroid) of the triangle
    double x_center = (a[0] + b[0] + c[0]) / 3;
    double y_center = (a[1] + b[1] + c[1]) / 3;
    double[] center = {x_center, y_center};

    
    //get the coordinates of the current point in the grid
    double[] currPoint = {voxelValue, gradMagnitude};

    // if the current point is not inside the triangle
    if (!isInside(a[0],a[1], b[0], b[1], c[0], c[1], currPoint[0], currPoint[1])){
        opacity = 0;
    }
    else{ // if the current point is  inside the triangle
          double y = currPoint[1];
          double x = currPoint[0];
          if (x == material_value){ // if voxel value is equal to the material value
            opacity = 1;
          }
          else if (x > material_value){ // if it is on the right part of the triangle
              
              // calculate the distance form the marerial value and the maximum value of this distance
              double max_dist =  calculateDistanceFromLineXaxis(b[0], b[1], c[0], c[1], y, material_value);
              double distance  = Math.abs(x - material_value);
        
    //        opacity is linearly decreasing as the distance from the material value increases
              opacity = 1 - (distance/max_dist);
          }
          else{ // if it is on the left part of the triangle
              
              // calculate the distance form the marerial value and the maximum value of this distance
              double max_dist =  calculateDistanceFromLineXaxis(a[0], a[1], c[0], c[1], y, material_value);
              double distance  = Math.abs(x - material_value);
              
    //        opacity is linearly decreasing as the distance from the material value increases
              opacity = 1 - (distance/max_dist);
          }
        }
    
    return opacity;
    }  


    // a utility function to calculate area of triangle  
    //formed by (x1, y1) (x2, y2) and (x3, y3)
    static double area(double x1, double y1, double x2, double y2, 
                                        double x3, double y3) 
    { 
       return Math.abs((x1*(y2-y3) + x2*(y3-y1)+ 
                                    x3*(y1-y2))/2.0); 
    } 
       
    // a function to check whether point P(x, y) lies 
    // inside the triangle formed by A(x1, y1), 
    // B(x2, y2) and C(x3, y3) */
    static boolean isInside(double x1, double y1, double x2, 
                double y2, double x3, double y3, double x, double y) 
    {    
       // calculate area of triangle ABC
        double A = area (x1, y1, x2, y2, x3, y3); 
       
        // calculate area of triangle PBC  
        double A1 = area (x, y, x2, y2, x3, y3); 
       
        // calculate area of triangle PAC
        double A2 = area (x1, y1, x, y, x3, y3); 
       
        // calculate area of triangle PAB */   
        double A3 = area (x1, y1, x2, y2, x, y); 
         
        // check if sum of A1, A2 and A3 is same as A
        return (A == A1 + A2 + A3); 
    }
    
    // a function to calculate the Euclidean distance between two points
    public static double calculateDistance(double[] A, double[] B)
    {
        double sum = 0;
        for (int i=0;i<A.length;i++){
            sum += Math.pow((A[i] - B[i]), 2.0); 
        }
    return Math.sqrt(sum);
    }
    
    // a function that calculates the distance  between vixel value and the 
    // maximum value the intensity can have inside the triangle for the given 
    // gradient of the voxel
    public static double calculateDistanceFromLineXaxis(double x1, double y1, double x2,double y2, double y, double material_value){
        double x = (y*(x2-x1) + x1*y2 - y1*x2)/ (y2-y1);
        return Math.abs(x - material_value);
    
    }

