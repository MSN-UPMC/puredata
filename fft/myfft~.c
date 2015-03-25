#include <stdlib.h>
#include "myfft_fft.c"
#include "myfft~.h"

#define PI 3.14159265358979323846
#define DSP_ADD_LENGTH 4
#define VECTOR_SIZE 4096
#define WAIT_UNTIL_PERFORM 4

/* Fill vector with begin values according to BlackMan
 * @param A vector to fill
 */
void init_with_blackman(t_sample *in){

  int i;
  for (i=0;i<VECTOR_SIZE;i++) {
    in[i]=(float)(0.42-0.5*(cos(2*PI*i/VECTOR_SIZE))+0.08*(cos(4*PI*i/VECTOR_SIZE)));
  }

}

/* Apply BlackMan to a vector
 * @param A vector to be applied
 * @param A vector with blackman values
 */
void apply_blackman(t_sample *in, t_sample *blackman){
  int i=0;
  for(i=0;i<VECTOR_SIZE;i++){
    in[i]=blackman[i]*in[i];
  }
}

/* Copy a signal into another
 * @param A vector model
 * @param A vector copied from the model
 * @param Size of vector 
 */
void copy_to_out_signal(t_sample *in, t_sample *out, int n){
  int i=0;
  for(i=0;i<n;i++){
    out[i]=in[i];
  }  
}

/* Shift n times on a signal
 * @param A vector to be shift
 * @param Shift number
 */
void shift_signal(t_sample *in, int n){
  int i=0;
  for(i=0;i<VECTOR_SIZE-n;i++){
    in[i]=in[i+n];
  }
}

/* Render a fft by using a circular buffer
 * @param incoming signal to be applied
 * @param blackman bitshuffle
 * @param blackman weighting
 */
void apply_buffered_fft(t_sample *in, int *bitshuffle, float *weighting){
  int i=0;  
  for(i=0;i<=VECTOR_SIZE;i+=1024){
    rdft((VECTOR_SIZE/WAIT_UNTIL_PERFORM),1,in+i,bitshuffle,weighting);
  }
}

/* Clear the buffer and add a signal into this one
 * @param buffer
 * @param signal input
 * @param size of signal
 */
void refresh_buffer_with_incoming_signal(t_sample *buffer, t_sample*in, int n){
  int i=0;
  for(i=0;i<n;i++) {
    buffer[VECTOR_SIZE-n+i]=*in++;
  }
}

/* Trigger function for Pure Data DSP */
t_int *myfft_tilde_perform(t_int *w){

  // Fetch objects
  t_myfft_tilde *x=(t_myfft_tilde *)(w[1]);
  t_sample *in =(t_sample *)(w[2]);
  t_sample *out =(t_sample *)(w[3]);  
  int n =(int)(w[4]);
  
  // refresh data on which work
  refresh_buffer_with_incoming_signal(x->circularbuffer,in,n);
  x->bitnumber+= VECTOR_SIZE;
  
  // If there is enough bits
  if(x->bitnumber>=(WAIT_UNTIL_PERFORM*VECTOR_SIZE)){
    x->bitnumber=0;
    x->complete=1;
    // Apply fft
    apply_blackman(x->circularbuffer,x->window);
    apply_buffered_fft(x->circularbuffer,x->bitshuffle,x->weighting);
  }
  // Modify output signal
  if(x->complete==1){
    copy_to_out_signal(x->circularbuffer,out,n);
  }
  // Remove n first bit from circular buffer
  shift_signal(x->circularbuffer,n);
  
  return (w+DSP_ADD_LENGTH+1);
}

/* Add myfft~ into DSP stack */
void myfft_tilde_dsp(t_myfft_tilde *x, t_signal **sp){
  dsp_add(myfft_tilde_perform, DSP_ADD_LENGTH, x,
    sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}

/* Free the memory of myfft~ */
void myfft_tilde_free(t_myfft_tilde *x){
  free(x->window);
  free(x->bitshuffle);
  free(x->weighting);
  free(x->circularbuffer);
  outlet_free(x->x_out);
}

/* Initialize myfft~ on PureData */
void *myfft_tilde_new(void){

  t_myfft_tilde *m = (t_myfft_tilde *)pd_new(myfft_tilde_class);
  m->x_out = outlet_new(&m->x_obj, &s_signal);
  
  // Init structure
  m->window = malloc(VECTOR_SIZE*(sizeof*m->window));
  m->bitshuffle = malloc(2*VECTOR_SIZE*(sizeof*m->bitshuffle));
  m->weighting = malloc(2*VECTOR_SIZE*(sizeof*m->weighting));
  m->circularbuffer = malloc(VECTOR_SIZE*(sizeof*m->circularbuffer));
  
  // Init default values
  m->complete=0;
  m->bitnumber=0;
  init_with_blackman(m->window);
  init_rdft(VECTOR_SIZE, m->bitshuffle, m->weighting);
  
  return (void *)m;
}

/* Setup myfft~ on PureData */
void myfft_tilde_setup(void){

  myfft_tilde_class = class_new(gensym("myfft~"),(t_newmethod)myfft_tilde_new,(t_method) myfft_tilde_free,sizeof(t_myfft_tilde),CLASS_DEFAULT,0);
  class_addmethod(myfft_tilde_class,(t_method)myfft_tilde_dsp,gensym("dsp"),0);
  CLASS_MAINSIGNALIN(myfft_tilde_class,t_myfft_tilde,f);
}
