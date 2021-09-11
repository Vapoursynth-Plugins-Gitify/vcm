#ifndef UnitSq2Quad_V_C_MOHAN
#define UnitSq2Quad_V_C_MOHAN
/*
Given 4 corners coordinates x0 y0,  x1 y1, x2 y2, x3 y3 
of the quadrilateral creates a 3 x 3 matrix for conversion 
from unit square to the quadrilateral.
 Also creates an inverse matrix inv for quad to 
 unit square conversion
sq is 3 X 3
quad is 4 X 2

  inv = 3 * 3

  returns -1 in case of error or noninvertible matrix
  otherwise returns 0
  */
int UnitSq2Quad(float sq[][3], float inv[][3], float quad[][2]);
//...............................................
int UnitSq2Quad(float sq[][3], float inv[][3], float quad[][2])
{
	
	float px = quad[0][0]- quad[1][0] + quad[2][0] - quad[3][0];

	float py = quad[0][1]- quad[1][1] + quad[2][1] - quad[3][1];

	if ( (px > -1e-10 && px < 1e-10) && (py > -1e-10 && py < 1e-10))
	{
		sq[0][0] = quad[1][0] - quad[0][0];

		sq[1][0] = quad[2][0] - quad[1][0];

		sq[2][0] = quad[0][0];

		sq[0][1] = quad[1][1] - quad[0][1];

		sq[1][1] = quad[2][1] - quad[1][1];

		sq[2][1] = quad[0][1];

		sq[0][2] = 0;

		sq[1][2] = 0;

		sq[2][2] = 1;

	}
	else
	{
		float dx1 = quad[1][0] - quad[2][0];
        float dx2 = quad[3][0] - quad[2][0];
        float dy1 = quad[1][1] - quad[2][1];
        float dy2 = quad[3][1] - quad[2][1];

		float det = dx1 * dy2 - dx2 * dy1;

		if ( det < -1e-10 || det > 1e-10)
		{
			sq[0][2] = (px * py - py * dx2) / det;

			sq[1][2] = (py * dx1 - px * dy1) / det;

			sq[2][2] = 1;

			sq[0][0] = quad[1][0] - quad[0][0] + sq[0][2]*quad[1][0];

			sq[1][0] = quad[3][0] - quad[0][0] + sq[1][2]*quad[3][0];

			sq[2][0] = quad[0][0];

			sq[0][1] = quad[1][1] - quad[0][1] + sq[0][2]*quad[1][1];

			sq[1][1] = quad[3][1] - quad[0][1] + sq[1][2]*quad[3][1];

			sq[2][1] = quad[0][1];
		}

		else
		{
			for(char i = 0; i < 3; i ++)

				for( char j = 0; j < 3; j++)

					sq[i][j] = 0;

			sq[0][0] = 1;
			sq[1][1] = 1;

			return -1;
		
		}
	}
			
		// invert the sq matrix into Inverted SQ matrix
		// This enables quad to unit square transform
/*       0 1 2
		|a b c| 0
	M = |d e f| 1
		|g h k| 2

	det(M) = a * ( e * k - f * h) - b * ( d * k - f * g) + c * ( d * h - e * g);
			= sq[0,0]* A = ([1,1] * [2,2] - [1,2] * [2,1])
			- sq[0,1]* B = ([1,0] * [2,2] - [1,2] * [2,0])
			+ sq[0,2]* C = ([1,0] * [2,1] - [1,1] * [2,0])

	det(M) should be non zero 

					|A D G |
	M- = 1 / det(M)	|B E H |
					|C F K |

	where A = ek-fh, B = fg-dk, C = dh-eg,	D = ch - bk, E = ak - cg, F = gb - ah
	G = bf - ce, H = cd - af, K = ae - bd
	A = inv[0,0] = [1,1] * [2,2] - [1,2] * [2,1]
	B = inv[1,0] = [1,0] * [2,2] - [1,2] * [2,0]
	C = inv[2,0] = [1,0] * [2,1] - [1,1] * [2,0]
	D = inv[0,1] = [0,2] * [2,1] - [0,1] * [2,2]
	E = inv[1,1] = [0,0] * [2,2] - [0,2] * [2,0]
	F = inv[2,1] = [2,0] * [0,1] - [0,0] * [2,1]
	G = inv[0,2] = [0,1] * [1,2] - [0,2] * [1,1]
	H = inv[1,2] = [0,2] * [1,0] - [0,0] * [1,2]
	K = inv[2,2] = [0,0] * [1,1] - [0,1] * [1,0]
*/

	float det = sq[0][0] * ( sq[1][1] * sq[2][2] - sq[1][2] * sq[2][1])
	  
		- sq[0][1] * ( sq[2][2] * sq[1][0] - sq[1][2] * sq[2][0])
	  
		+ sq[0][2] * ( sq[1][0] * sq[2][1] - sq[1][1] * sq[2][0]);
	// non zero det
	if ( det < -1e-10 || det > 1e-10)
	{
	  inv[0][0] = ( sq[1][1] * sq[2][2] - sq[1][2] * sq[2][1]) / det;
	  inv[1][0] = ( sq[1][2] * sq[2][0] - sq[1][0] * sq[2][2]) / det;
	  inv[2][0] = ( sq[1][0] * sq[2][1] - sq[1][1] * sq[2][0]) / det;
	  inv[0][1] = ( sq[0][2] * sq[2][1] - sq[0][1] * sq[2][2]) / det;
	  inv[1][1] = ( sq[0][0] * sq[2][2] - sq[0][2] * sq[2][0]) / det;
	  inv[2][1] = ( sq[2][0] * sq[0][1] - sq[0][0] * sq[2][1]) / det;
	  inv[0][2] = ( sq[0][1] * sq[1][2] - sq[0][2] * sq[1][1]) / det;
	  inv[1][2] = ( sq[0][2] * sq[1][0] - sq[0][0] * sq[1][2]) / det;
	  inv[2][2] = ( sq[0][0] * sq[1][1] - sq[0][1] * sq[1][0]) / det;
	}

	else
		return -1;

	return 0;
 }

	
#endif




