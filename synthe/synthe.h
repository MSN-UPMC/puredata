#ifndef _synthe_h_
#define _synthe_h_

#include "m_pd.h"

static t_class *synthe_class;

typedef struct _synthe {
  t_object x_obj;

  t_sample f;

  t_inlet *modulatrice;
  t_inlet *messages;
  t_outlet *x_out;

  // BlackMan usefull structure
  int *bitshuffle;
  float *weighting;
  float *window;
  
  // Synthe usefull structure
  t_int bypass;
  t_int autonorm;
  t_float shapeWidth;
} t_synthe;

/* Trigger function for Pure Data DSP 
 *  @param  Vector which contains DSP add informations 
 *  @return Size of Vector and add one 
 */
t_int *synthe_perform(t_int *w);

/* Add myfft~ into DSP stack 
 *  @param myfft~ object
 *  @param input signal
 */
void synthe_dsp(t_synthe *x, t_signal **sp);

/* Free the memory of myfft~ 
 *  @param myfft~ object 
 */
void synthe_free(t_synthe *x);

/* Initialize synthe on PureData 
   @param (autonorm value)? (bypass value)?
 */
void *synthe_new(int argc, t_atom *argv);

/* Setup synthe on PureData */
void synthe_setup(void);

/* Message inlet to manage autonorm and bypass 
 * @param autonorm value
 * @param bypass value
*/
void synthe_messages(t_synthe *x, t_floatarg n, t_floatarg p);

#endif
