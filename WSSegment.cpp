/********************************************************************************
WSSegment segments an image by Allen Soille algorithm of water shed segmentation.
Basin numbers and watersheds are determined by this function and marked in tag 
buffer basin number and in  ws buffer as zero for watershed..

  All buffers must be of size height * width.

  The input is pointers to image values sorted in ascending order of image values
  to which they point, pointers to buffer where image resides,
  tag in which basin numbers will be computed, dist for indicating distance from minima,
  and ws where watershed will be recorded as one and rest as 0. image width, image height,

Author V.C.Mohan
Oct 22, 2006, mar 24, 2009, oct 2014(added template)
*********************************************************************************/  
#ifndef VCM_WSSEGMENT_CPP
#define VCM_WSSEGMENT_CPP   
#include <queue>
template <typename finc>
void WSSegment( const finc **pixelsort, const finc *buf, int * tag,
		  int * dist, char * ws, int ht, int wd, bool connect4 = true);
template <typename finc>
void WSSegment( const finc **pixelsort, const finc *buf, int * tag,
		  int * dist, char * ws, int ht, int wd, bool connect4)
{


				std::queue <const finc *> fifoq;

					// initialize 
				int init=-1, mask=-2,wshed=0;
				int curtag=0, curdist=0;
				int psize = ht * wd;
				
				// prepare data for sorting

				for(int h=0; h < psize; h++)
					
				{
					
					tag[h] = init;
					dist[h]= 0;
					ws[h] = init;
				}

						// minimum and max grey values
				int gmin = *(pixelsort[0]) ;
			
				int gmax = *(pixelsort[psize-1]);
			
				int gindex = 0, pindex = 0;		// these will hopefully reduce search of pixelsort buffer for a g value
					// start flooding

				for( int g = gmin; g <= gmax; g++)	// geodesic SKIZ of level g-1 inside level g
				{
					// get pixels having this g value
				for(int j = gindex; j < psize; j++)
				{
					if(*pixelsort[j]  < g )
						continue;		// we already completed this value

					if(*pixelsort[j] > g )
						
						break;		// we reached next value

						// same value pixel. get its original coordinates
				//	int h = (pixelsort[j] - buf);
					

					tag[pixelsort[j] - buf] = mask;		//  mask all pixels with this g value. Do I need to mark all prior to rest of code?
												// yes. done it by seperating loops
				}

				for(int j = gindex; j < psize; j++)
				{
					if(*pixelsort[j] < g )
						continue;
					if(*pixelsort[j] > g )
					{
						gindex=j;
						break;
					}
						// same value pixel. get its original coordinates
					int h = (pixelsort[j] - buf) / wd;

					int w = (pixelsort[j] - buf) % wd;

					int hwd = h * wd;

					if(h > 0 && h < ht - 1 && w > 0 && w < wd - 1)	// if neighbour is in frame
					{
							// check whether 4 close neighbour tag is >0 or is wshed
						if(	tag[hwd + wd + w] >0	|| ws[hwd + wd + w] == wshed ||
							tag[hwd + w -1] >0	|| ws[hwd + w - 1] == wshed ||
							tag[hwd + w +1] >0	|| ws[hwd + w + 1] == wshed ||
							tag[hwd - wd + w] >0	|| ws[hwd - wd + w] == wshed	)

						{
							dist[hwd + w]=1;
								// initialize queue with neighbours at level g of current basins or water sheds

							fifoq.push(pixelsort[j]);
						}
						 else if(!connect4)		// 8 connect option. so check rest four distant neighbours (at corners) also
						{
							if(	tag[hwd + wd+w-1] > 0 || ws[hwd + wd+w-1] == wshed ||
								tag[hwd - wd+w-1] > 0 || ws[hwd - wd+w-1] == wshed ||
								tag[hwd + wd+w+1] > 0 || ws[hwd + wd+w+1] == wshed ||
								tag[hwd - wd+w+1] > 0 || ws[hwd - wd+w+1] == wshed 	)
							{
								dist[hwd + w] = 1;

									// initialize queue with neighbours at level g of current basins or water sheds
								fifoq.push(pixelsort[j]); 
							}
						}

							
					}
				}


					// up date curdist
				curdist=1;

				fifoq.push(NULL);		// a fictitous value used as a control marker

					// indefinite loop. breaks when queue is empty

				while(true)		// extend basins
				{
					const finc * p = fifoq.front();
					fifoq.pop();
									// if not a real value
					if(p == NULL)
					{
						if(fifoq.empty())
							break;
						else
						{
							curdist++;
							fifoq.push(NULL);	// another marker
							p = fifoq.front();
							fifoq.pop();
						}
					}

						// get back x, y coordinates and its value
					int h = (p - buf) / wd;
					int w = (p - buf) % wd;


						// check all connected neighbours of point p on dist and tag frame

					for(int s = -1; s < 2; s++)
						for(int t = -1; t < 2; t++)
						{
							// check for 4 connectivity
							if (connect4)
							{
								if(s == t || s==-t)		// not close neighbours
									continue;
							}
							// check for 8 connectivity 
							else
								if(s == 0 && t == 0)		// same point
									continue;

								
							if(h + s < 0 || h + s > ht - 1 || w + t < 0 || w + t > wd - 1)
								continue;		// not in frame. so skip this point

							int hwd = h * wd;

							int hswd = (h + s) * wd;

							if(dist[hswd + (w + t)] < curdist 
								&& (tag[hswd + (w + t)] > 0 
								|| ws[hswd + (w + t)] == wshed) )
							{
									// this neighbour point belongs to existing basin or watershed
								if(tag[hswd + (w + t)] > 0)
								{
									//  yes this neighbour belongs to an existing basin

									if(tag[hwd + w] == mask || ws[hwd + w] == wshed)
											// put same basin tag on p as that of neighbour 
										tag[hwd + w] = tag[hswd + (w+t)];

									else if( tag[hwd + w] != tag[hswd + (w+t)])
										// this neighbour belongs to a diff basin. so our point is water shed
										ws[hwd + w] = wshed;
								}

								else
									if( tag[hwd + w] == mask)
										ws[hwd + w] = wshed;
							}


							else if( tag[hswd + (w + t)] == mask 
									&& dist[hswd + (w + t)] == 0)
							{
										// this neighbour is a plataue pixel
								dist[hswd + (w + t)] = curdist + 1;

								fifoq.push( buf + hswd + (w + t) );
							}
						}
					
				}	// while(true)


					// once again  find new minima at this same level g
					// get pixels having this g value
					// pindex was originally 0.
				for(int j = pindex; j < psize; j++)
				{
					if(*pixelsort[j] < g )
						continue;
					if(*pixelsort[j] > g )
					{
						pindex=j;
						break;
					}
						// same value pixel
					int h = (pixelsort[j] - buf) / wd;

					int w = (pixelsort[j] - buf) % wd;

					dist[pixelsort[j] - buf] = 0;	// reset to zero

					if( h >= 0 && h < ht && w >= 0 && w < wd)	// is this reqd? why h>0 and not h>=0 ?
					{
						if(tag[h * wd + w] == mask)		// point is inside new minimum
							// check whether 4 neighbour tag is >0 or is wshed
						{
							curtag++;		// create new label 

							fifoq.push (pixelsort[j]);

							tag[h * wd + w] = curtag;

							while( ! fifoq.empty())
							{
								const finc * q = fifoq.front();
								fifoq.pop();

								int qh = (q - buf) / wd;

								int qw = (q - buf) % wd;


								for(int qs = -1; qs < 2; qs++)
									for(int qt = -1; qt < 2; qt++)
									{
										if(connect4)
										{
											if( qs == qt || qs == -qt)
												continue;
										}
										else
											if( qs == 0 && qt == 0)
												continue;

										if( qs + qh < 0 || qs + qh >= ht || qw + qt < 0 || qw + qt >= wd)
											continue;
												// check the neighbours of this point q
										if(tag[(qs+qh) * wd + (qw+qt)] == mask)
										{
											fifoq.push( buf + (qs+qh) * wd +  (qw+qt) );
											tag[(qs+qh) * wd + (qw+qt)] = curtag;
										}
									}
							}	// while( ! fifoq.empty())

						}	// if(tag[h * wd + w] == mask)

					}	//if( h > 0 && h < bht && w > 0 && w < bwd)

				}	// for(j=pindex;j<psize;j++)

			}	// for( int g = gmin; g <= gmax; g++)

		}
#endif		
//--------------------------------------------------------------------------------------------------		
