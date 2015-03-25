#ifndef _fft_h_
#define _fft_h_

#include "m_pd.h"

static t_class *myfft_tilde_class;

typedef struct _myfft_tilde {
  t_object x_obj;

  t_sample f_pan;
  t_sample f;
  t_sample *circularbuffer;
  
  t_outlet *x_out;

  // BlackMan usefull structure
  int *bitshuffle;
  float *weighting;
  float *window;
  int bitnumber;
  int complete;

} t_myfft_tilde;

/* Trigger function for Pure Data DSP 
 *  @param  Vector which contains DSP add informations 
 *  @return Size of Vector and add one 
 */
t_int *fft_tilde_perform(t_int *w);

/* Add myfft~ into DSP stack 
 *  @param myfft~ object
 *  @param input signal
 */
void myfft_tilde_dsp(t_myfft_tilde *x, t_signal **sp);


/* Free the memory of myfft~ 
 *  @param myfft~ object 
 */
void myfft_tilde_free(t_myfft_tilde *x);

/* Initialize myfft~ on PureData */
void *myfft_tilde_new(void);

/* Setup myfft~ on PureData */
void myfft_tilde_setup(void);

#endif
