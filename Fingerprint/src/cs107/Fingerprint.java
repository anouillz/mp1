package cs107;

import java.util.ArrayList;

import java.util.List;

/**
 * Provides tools to compare fingerprint.
 */
public class Fingerprint {

	/**
	 * The number of pixels to consider in each direction when doing the linear
	 * regression to compute the orientation.
	 */
	public static final int ORIENTATION_DISTANCE = 16;

	/**
	 * The maximum distance between two minutiae to be considered matching.
	 */
	public static final int DISTANCE_THRESHOLD = 5;

	/**
	 * The number of matching minutiae needed for two fingerprints to be considered
	 * identical.
	 */
	public static final int FOUND_THRESHOLD = 20;

	/**
	 * The distance between two angle to be considered identical.
	 */
	public static final int ORIENTATION_THRESHOLD = 20;

	/**
	 * The offset in each direction for the rotation to test when doing the
	 * matching.
	 */
	public static final int MATCH_ANGLE_OFFSET = 2;

	/**
	 * Returns an array containing the value of the 8 neighbours of the pixel at
	 * coordinates <code>(row, col)</code>.
	 * <p>
	 * The pixels are returned such that their indices corresponds to the following
	 * diagram:<br>
	 * ------------- <br>
	 * | 7 | 0 | 1 | <br>
	 * ------------- <br>
	 * | 6 | _ | 2 | <br>
	 * ------------- <br>
	 * | 5 | 4 | 3 | <br>
	 * ------------- <br>
	 * <p>
	 * If a neighbours is out of bounds of the image, it is considered white.
	 * <p>
	 * If the <code>row</code> or the <code>col</code> is out of bounds of the
	 * image, the returned value should be <code>null</code>.
	 *
	 * @param image array containing each pixel's boolean value.
	 * @param row   the row of the pixel of interest, must be between
	 *              <code>0</code>(included) and
	 *              <code>image.length</code>(excluded).
	 * @param col   the column of the pixel of interest, must be between
	 *              <code>0</code>(included) and
	 *              <code>image[row].length</code>(excluded).
	 * @return An array containing each neighbours' value.
	 */
	public static boolean[] getNeighbours(boolean[][] image, int row, int col) {
		assert (image != null); // special case that is not expected (the image is supposed to have been checked
		// earlier)

		boolean[] retour = new boolean[8];

		if ( (row < 0) || (row>= image.length)||(col<0)||(col>=image[0].length) ) {

			return null ;
		}

		if(col == 0){
			retour[5] = false;
			retour[6] = false;
			retour[7] = false;
		}else{
			retour[6] = image[row][col-1];
		}
		if(col == image[0].length - 1){
			retour[1] = false;
			retour[2] = false;
			retour[3] = false;
		}else{
			retour[2] = image[row][col+1];
		}
		if(row == 0){
			retour[7] = false;
			retour[0] = false;
			retour[1] = false;
		}else{
			retour[0] = image[row-1][col];
		}
		if(row == image.length - 1){
			retour[5] = false;
			retour[4] = false;
			retour[3] = false;
		}else{
			retour[4] = image[row+1][col];
		}

		if(col != 0 && row != 0){
			retour[7] = image[row-1][col-1];
		}
		if(row != 0 && col != image[0].length - 1){
			retour[1] = image[row-1][col+1];
		}
		if(col != image[0].length - 1 && row != image.length - 1){
			retour[3] = image[row+1][col+1];
		}
		if(col != 0 && row != image.length - 1){
			retour[5] = image[row+1][col-1];
		}

		return retour;
	}

	/**
	 * Computes the number of black (<code>true</code>) pixels among the neighbours
	 * of a pixel.
	 *
	 * @param neighbours array containing each pixel value. The array must respect
	 *                   the convention described in
	 *                   {@link #getNeighbours(boolean[][], int, int)}.
	 * @return the number of black neighbours.
	 */
	public static int blackNeighbours(boolean[] neighbours) {

		int compteurBlackNeighbours = 0 ;



		for (int i = 0 ; i < neighbours.length ; ++i ) {

			if (neighbours[i] == true) {

				++compteurBlackNeighbours ;
			}
		}

		return compteurBlackNeighbours ;
	}

	/**
	 * Computes the number of white to black transitions among the neighbours of
	 * pixel.
	 *
	 * @param neighbours array containing each pixel value. The array must respect
	 *                   the convention described in
	 *                   {@link #getNeighbours(boolean[][], int, int)}.
	 * @return the number of white to black transitions.
	 */
	public static int transitions(boolean[] neighbours) {

		int foisDeTransition = 0 ;


		if ((neighbours[7] == false)&&(neighbours[0] == true)) {

			++foisDeTransition ;
		}


		for (int i = 1 ; i < neighbours.length ; ++i ) {


			if ((neighbours[i] == true)&&(neighbours [i-1] == false)) {

				++foisDeTransition ;

			} 

		}

		return foisDeTransition ;
	}

	/**
	 * Returns <code>true</code> if the images are identical and false otherwise.
	 *
	 * @param image1 array containing each pixel's boolean value.
	 * @param image2 array containing each pixel's boolean value.
	 * @return <code>True</code> if they are identical, <code>false</code>
	 *         otherwise.
	 */
	public static boolean identical(boolean[][] image1, boolean[][] image2) {


		int row;
		int col;



		if ((image1.length != image2.length)||(image1[0].length!=image2[0].length)) {

			return false ;
		}


		for(row = 0; row < image1.length; row++) {
			for (col = 0; col < image1[0].length; col++) {

				if (image1[row][col] != image2[row][col]) {

					return false ;
				}
			}
		}

		return true;
	}


	/**
	 * Internal method used by {@link #thin(boolean[][])}.
	 *
	 * @param image array containing each pixel's boolean value.
	 * @param step  the step to apply, Step 0 or Step 1.
	 * @return A new array containing each pixel's value after the step.
	 */
	public static boolean[][] thinningStep(boolean[][] image, int step) {

		int row;
		int col;
		int imageBlackNeighbours;
		int imageTransitions;




		boolean[] imageNeighbours = new boolean [8];

		boolean[][] newImage = new boolean [image.length][image[0].length];




		if ((step!=0)&&(step!=1)) {

			return null ;
		}

		for(row = 0; row < image.length; row++) {

			for (col = 0; col < image[0].length; col++) {

				newImage[row][col] = image[row][col] ;
				imageNeighbours = getNeighbours(image,row,col);
				imageBlackNeighbours = blackNeighbours(imageNeighbours);
				imageTransitions = transitions(imageNeighbours);

				if (image[row][col] == true) {
					if (imageNeighbours != null) {
						if ((imageBlackNeighbours >= 2) && (imageBlackNeighbours <= 6)) {
							if (imageTransitions == 1) {

								if (step == 0) {
									if  ((imageNeighbours[0] == false) || (imageNeighbours[2] == false) || (imageNeighbours[4] == false)) {
										if ((imageNeighbours[2] == false) || (imageNeighbours[4] == false) || (imageNeighbours[6] == false)) {

											newImage[row][col] = false;



										}

									}

								}else {
									if  ((imageNeighbours[0] == false) || (imageNeighbours[2] == false) || (imageNeighbours[6] == false)) {
										if ((imageNeighbours[0] == false) || (imageNeighbours[4] == false) || (imageNeighbours[6] == false)) {

											newImage[row][col] = false;



										}

									}



								}
							}
						}


					}
				}
			}

		}


		return newImage;
	}

	/**
	 * Compute the skeleton of a boolean image.
	 *
	 * @param image array containing each pixel's boolean value.
	 * @return array containing the boolean value of each pixel of the image after
	 *         applying the thinning algorithm.
	 */
	public static boolean[][] thin(boolean[][] image) {

		boolean[][] newImage = new boolean [image.length][image[0].length];
		boolean[][] oldImage = new boolean [image.length][image[0].length];

		do {

			oldImage = newImage ;			
			newImage = thinningStep(image , 0) ;
			newImage = thinningStep(newImage , 1) ;

			

		}while (!identical(oldImage,newImage)) ;

		return newImage;
	}

	/**
	 * Computes all pixels that are connected to the pixel at coordinate
	 * <code>(row, col)</code> and within the given distance of the pixel.
	 *
	 * @param image    array containing each pixel's boolean value.
	 * @param row      the first coordinate of the pixel of interest.
	 * @param col      the second coordinate of the pixel of interest.
	 * @param distance the maximum distance at which a pixel is considered.
	 * @return An array where <code>true</code> means that the pixel is within
	 *         <code>distance</code> and connected to the pixel at
	 *         <code>(row, col)</code>.
	 */
	public static boolean[][] connectedPixels(boolean[][] image, int row, int col, int distance) {

		int connectedPixelsBlackNeighbours ;
		boolean[]connectedPixelsNeighbours ;

		boolean[][] connectedPixels = new boolean [image.length][image[0].length];
		connectedPixels[row][col]=true;



		int row1 ;
		int col1 ;


		if( !(row<0 ||col <0 || row>=image.length || col >= image[0].length)) {
			boolean test = true  ;
			while (test) {

				test=false;

				for(row1 = 0; row1 < image.length; row1++) {

					for (col1 = 0; col1 < image[0].length; col1++) {


						connectedPixelsNeighbours = getNeighbours(connectedPixels,row1,col1);
						connectedPixelsBlackNeighbours = blackNeighbours(connectedPixelsNeighbours);	

						if (image[row1][col1]==true) {
							//pixel est noir?


							if (connectedPixels[row1][col1] != true) {
								//pixel est déjà dans le tableau?


								if (connectedPixelsBlackNeighbours>0){
									//autres voisins noirs?

									if((row1<=row+distance)&&(col1<=col+distance)&&(row1>=row-distance)&&(col1>=col-distance)) {

										test = true ;



										connectedPixels[row1][col1]=true ;



									}

								}

							}



						}


					}


				}
			}
		}


		return connectedPixels;
	}

	/**
	 * Computes the slope of a minutia using linear regression.
	 *
	 * @param connectedPixels the result of
	 *                        {@link #connectedPixels(boolean[][], int, int, int)}.
	 * @param row             the row of the minutia.
	 * @param col             the col of the minutia.
	 * @return the slope.
	 */
	public static double computeSlope(boolean[][] connectedPixels, int row, int col) {
		// row of the green pixel
		int row1;

		// col of the green pixel
		int col1; 

		// variable that represents the value of the slope
		double a = 0.0;

		double sum = 0;

		double sumy = 0;

		double sumx = 0;

		for (row1 = 0; row1 < connectedPixels.length; row1++) {
			for(col1 = 0; col1 < connectedPixels[0].length; col1++) {

				// x and y the coordinates of the black pixel
				int x = col1 - col;
				int y = row - row1;


				// if the pixel is black 
				if (connectedPixels[row1][col1] == true) {


					// gives the sum of x*y

					sum = sum + (x*y);

					// gives the sum of x^2

					sumx = sumx + (x*x);

					//gives the sum of y^2

					sumy = sumy + (y*y);

				}
			}
		}

		// different value for a depending on the value of sum x^2 and sum y^2
		if (sumx >= sumy) {

			if (sumx == 0) {

				return Double.POSITIVE_INFINITY;

			}

			a = (sum / sumx);

			return a;

		} else {



			a = (sumy / sum);

			return a;



		}


	}


	/**
	 * Computes the orientation of a minutia in radians.
	 * 
	 * @param connectedPixels the result of
	 *                        {@link #connectedPixels(boolean[][], int, int, int)}.
	 * @param row             the row of the minutia.
	 * @param col             the col of the minutia.
	 * @param slope           the slope as returned by
	 *                        {@link #computeSlope(boolean[][], int, int)}.
	 * @return the orientation of the minutia in radians.
	 */
	public static double computeAngle(boolean[][] connectedPixels, int row, int col, double slope) {



		double angle = 0 ;


		double p = Math.PI;

		//calculates number of pixels that are above or under the line
		int nb_dessus = 0;
		int nb_dessous = 0;


		for (int row1 = 0; row1 < connectedPixels.length ; row1++) {
			for(int col1 = 0; col1 < connectedPixels[0].length ; col1++) {

				// x and y the coordinates of the black pixel
				int x = col1 - col;
				int y = row - row1;

				if (connectedPixels[row1][col1]==true) {

					// equation of the line
					if (y >= (-x/slope)) {

						nb_dessus++;

					}else {

						nb_dessous++;

					}
				}
			}
		}



		// value of the angle if the slope = double.positive_infinity
		if ( slope == Double.POSITIVE_INFINITY) {


			if (nb_dessous >= nb_dessus) {

				angle = -p/2;

			} else {


				angle = p/2;

			}

			return angle;

		}

		angle = Math.atan(slope);

		if (angle >= 0 && nb_dessous >= nb_dessus) {

			angle = angle + p ;

			return angle;
		}


		if ( angle < 0 && nb_dessous < nb_dessus) {

			angle = angle + p ;

			return angle;





		}

		return angle;
	}

	/**
	 * Computes the orientation of the minutia that the coordinate <code>(row,
	 * col)</code>.
	 *
	 * @param image    array containing each pixel's boolean value.
	 * @param row      the first coordinate of the pixel of interest.
	 * @param col      the second coordinate of the pixel of interest.
	 * @param distance the distance to be considered in each direction to compute
	 *                 the orientation.
	 * @return The orientation in degrees.
	 */
	public static int computeOrientation(boolean[][] image, int row, int col, int distance) {


		boolean [][] imageConnected = new boolean[image.length][image[0].length];
		imageConnected = connectedPixels(image, row, col, distance);

		double slope = computeSlope(imageConnected, row, col);

		//angle in radiant
		double angle = computeAngle(imageConnected, row, col, slope);

		//angle but in degrees
		double ang = Math.toDegrees(angle);


		// o = orientation, the angle rounded to the closest integer
		int o = (int) Math.round(ang);


		// if orientation is strictly negative, we add 360
		if (o < 0) {

			o = o + 360;

		} 

		return o;

	}

	/**
	 * Extracts the minutiae from a thinned image.
	 *
	 * @param image array containing each pixel's boolean value.
	 * @return The list of all minutiae. A minutia is represented by an array where
	 *         the first element is the row, the second is column, and the third is
	 *         the angle in degrees.
	 * @see #thin(boolean[][])
	 */
	public static List<int[]> extract(boolean[][] image) {

		// row, colomn and orientation of the minutia
		int row;
		int col;
		int o;



		boolean [] imageNeighbours = new boolean [8];


		List<int[]> listOfM = new ArrayList <int[]>();	  


		for(row = 1; row < image.length - 1; row++) {
			for (col = 1; col < image[0].length - 1; col++) {




				// if the pixel is black and its transitions are equal to 1 or 3, we add the minutia to the list
				if (image[row][col] == true) {



					imageNeighbours = getNeighbours(image, row, col);

					// the number of transitions of the minutia
					int imageTransitions = transitions(imageNeighbours);



					if ((imageTransitions == 1)||(imageTransitions == 3)) {

						o = computeOrientation(image, row, col, ORIENTATION_DISTANCE);

						// array of integers of size 3, first element is row, second is col and third is orientation in degrees
						int[] minutia = new int[3];

						minutia[0] = row;
						minutia[1] = col;
						minutia[2] = o;

						listOfM.add(minutia);



					} 

				}

			}


		}

		return listOfM;

	}

	/**
	 * Applies the specified rotation to the minutia.
	 *
	 * @param minutia   the original minutia.
	 * @param centerRow the row of the center of rotation.
	 * @param centerCol the col of the center of rotation.
	 * @param rotation  the rotation in degrees.
	 * @return the minutia rotated around the given center.
	 */
	public static int[] applyRotation(int[] minutia, int centerRow, int centerCol, int rotation) {


		int row = minutia[0];
		int col = minutia[1];
		int o = minutia[2];

		int x = col - centerCol;
		int y = centerRow - row;

		double newRotation = Math.toRadians(rotation);

		double cosR = Math.cos(newRotation);
		double sinR = Math.sin(newRotation);

		double newX = x * cosR - y * sinR;
		double newY = x * sinR + y * cosR;

		double newRow = centerRow - newY;
		double newCol = newX + centerCol;
		int newOrientation = (o + rotation) % 360;

		int [] applyRotation = new int [3];

		int newR = (int) Math.round(newRow);
		int newC = (int) Math.round(newCol);


		applyRotation[0] = newR;
		applyRotation[1] = newC;
		applyRotation[2] = newOrientation;





		return (applyRotation);
	}

	/**
	 * Applies the specified translation to the minutia.
	 *
	 * @param minutia        the original minutia.
	 * @param rowTranslation the translation along the rows.
	 * @param colTranslation the translation along the columns.
	 * @return the translated minutia.
	 */
	public static int[] applyTranslation(int[] minutia, int rowTranslation, int colTranslation) {


		int row = minutia[0];
		int col = minutia[1];
		int o = minutia[2];

		double newRow = row - rowTranslation;
		double newCol = col - colTranslation;
		double newOrientation = o;

		int [] applyTranslation = new int [3];

		int newR = (int) Math.round(newRow);
		int newC = (int) Math.round(newCol);
		int newO = (int) Math.round(newOrientation);

		applyTranslation[0] = newR;
		applyTranslation[1] = newC;
		applyTranslation[2] = newO;



		return applyTranslation;
	} 

	/**
	 * Computes the row, column, and angle after applying a transformation
	 * (translation and rotation).
	 *
	 * @param minutia        the original minutia.
	 * @param centerCol      the column around which the point is rotated.
	 * @param centerRow      the row around which the point is rotated.
	 * @param rowTranslation the vertical translation.
	 * @param colTranslation the horizontal translation.
	 * @param rotation       the rotation.
	 * @return the transformed minutia.
	 */
	public static int[] applyTransformation(int[] minutia, int centerRow, int centerCol, int rowTranslation,
			int colTranslation, int rotation) {
		int[] minutiaRotation = new int[3];
		int[] minutiaTranslation = new int[3];

		minutiaRotation = applyRotation(minutia, centerRow, centerCol, rotation);
		minutiaTranslation = applyTranslation(minutiaRotation, rowTranslation, colTranslation);

		return minutiaTranslation;
	}

	/**
	 * Computes the row, column, and angle after applying a transformation
	 * (translation and rotation) for each minutia in the given list.
	 *
	 * @param minutiae       the list of minutiae.
	 * @param centerCol      the column around which the point is rotated.
	 * @param centerRow      the row around which the point is rotated.
	 * @param rowTranslation the vertical translation.
	 * @param colTranslation the horizontal translation.
	 * @param rotation       the rotation.
	 * @return the list of transformed minutiae.
	 */
	public static List<int[]> applyTransformation(List<int[]> minutiae, int centerRow, int centerCol, int rowTranslation,
			int colTranslation, int rotation) {

		ArrayList <int[]> list = new ArrayList<int[]>();


		for (int i = 0; i < minutiae.size() - 1; i++) {

			list.set(i, applyTransformation(minutiae.get(i), centerRow, centerCol, rowTranslation, colTranslation, rotation));  

		}

		return list;
	}
	/**
	 * Counts the number of overlapping minutiae.
	 *
	 * @param minutiae1      the first set of minutiae.
	 * @param minutiae2      the second set of minutiae.
	 * @param maxDistance    the maximum distance between two minutiae to consider
	 *                       them as overlapping.
	 * @param maxOrientation the maximum difference of orientation between two
	 *                       minutiae to consider them as overlapping.
	 * @return the number of overlapping minutiae.
	 */
	public static int matchingMinutiaeCount(List<int[]> minutiae1, List<int[]> minutiae2, int maxDistance,
			int maxOrientation) {

		double deltaR;
		double deltaC;

		int count = 0;

		for(int i = 0; i < minutiae1.size(); i++) {
			for(int j = 0; j< minutiae2.size(); j++) {

				deltaR = minutiae1.get(i) [0] - minutiae2.get(j) [0];
				deltaC = minutiae1.get(i) [1] - minutiae2.get(j) [1];

				double distance = Math.sqrt((deltaR * deltaR) + (deltaC * deltaC));

				int diffO= Math.abs(minutiae1.get(i)[2] - minutiae2.get(j)[2]);

				if(distance <= maxDistance) {
					if (diffO <= maxOrientation) {

						count ++;
						
					}
				}
			}
		}

		return count;
	}

	/**
	 * Compares the minutiae from two fingerprints.
	 *
	 * @param minutiae1 the list of minutiae of the first fingerprint.
	 * @param minutiae2 the list of minutiae of the second fingerprint.
	 * @return Returns <code>true</code> if they match and <code>false</code>
	 *         otherwise.
	 */
	public static boolean match(List<int[]> minutiae1, List<int[]> minutiae2) {

		boolean value = false;

		ArrayList <int[]> newMinutiae2 = new ArrayList<int[]>();

		
		for (int i = 0; i < Math.min(minutiae1.size(), minutiae2.size()); i++) {
			

			int rowTrans = minutiae2.get(i)[0] - minutiae1.get(i)[0];
			int colTrans = minutiae2.get(i)[1] - minutiae1.get(i)[1];
			int rotation = minutiae2.get(i)[2] - minutiae1.get(i)[2];
			

			for (int x = rotation - MATCH_ANGLE_OFFSET; x <= rotation + MATCH_ANGLE_OFFSET; x++) {

				int [] transformedMinutiae = applyTransformation(minutiae2.get(i), minutiae1.get(i)[0], minutiae1.get(i)[1], rowTrans, colTrans, x);

				newMinutiae2.add(transformedMinutiae);
			
				
			}
		}

		int count = matchingMinutiaeCount(minutiae1, newMinutiae2, DISTANCE_THRESHOLD, ORIENTATION_THRESHOLD);

		if (count >= FOUND_THRESHOLD) {

			value = true;
		}


		return value;
	}
}


