// get interpolated value and store in appropriate format

#ifndef REFORMHELPER_CPP_V_C_MOHAN
#define REFORMHELPER_CPP_V_C_MOHAN
//--------------------------------------------------------------------------------------------------
template <typename finc>
void PutClampedInterpolatedValue( finc * dstp, const finc * srcp, int pitch, int plane,
								int h, int w, int sx, int sy,int qx, int qy, int span,float * lbuf, finc min, finc max);

template <typename finc>
void PutNearestPointvalue( finc * dstp, const finc * srcp, int pitch, int kb,
								int h, int w, int sx, int sy, finc min, finc max);
//---------------------------------------------------------------------------------------------------
template <typename finc>
void PutClampedInterpolatedValue( finc * dstp, const finc * srcp, int pitch, int plane,
								int h, int w, int sx, int sy,int qx, int qy, int span,float * lbuf, finc min, finc max)
{
							
	
		dstp[(h )  * pitch + w ] 
			= clamp( LaQuantile (srcp + sy * pitch + sx, pitch, span, qx, qy, lbuf), min, max);
	
}
//------------------------------------------------------------------------------------------------------
template <typename finc>
void PutNearestPointvalue( finc * dstp, const finc * srcp, int pitch, int kb,
								int h, int w, int sx, int sy)
{	
							
	
		dstp[h * pitch + w * kb ] 
			= srcp [ sy * pitch + sx * kb];
	
}

#endif