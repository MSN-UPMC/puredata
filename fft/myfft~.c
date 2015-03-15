#include <stdlib.h>
#include "myfft_fft.c"
#include "myfft~.h"

#define PI 3.14159265358979323846
#define DSP_ADD_LENGTH 4
#define VECTOR_SIZE 4096
#define WAIT_UNTIL_PERFORM 4


void init_with_blackman(t_sample *in){

  int i;
  for (i=0;i<VECTOR_SIZE;i++) {
    in[i]=(float)(0.42-0.5*(cos(2*PI*i/VECTOR_SIZE))+0.08*(cos(4*PI*i/VECTOR_SIZE)));
  }

}

void apply_blackman(t_sample *in, t_sample *blackman){
  int i=0;
  for(i=0;i<VECTOR_SIZE;i++){
    in[i]=blackman[i]*in[i];
  }
}

void copy_to_out_signal(t_sample *in, t_sample *out, int n){
  int i=0;
  for(i=0;i<n;i++){
    out[i]=in[i];
  }  
}

void shift_signal(t_sample *in, int n){
  int i=0;
  for(i=0;i<VECTOR_SIZE-n;i++){
    in[i]=in[i+n];
  }
}

void apply_buffered_fft(t_sample *in, int *bitshuffle, float *weighting){
  int i=0;  
  for(i=0;i<=VECTOR_SIZE;i+=1024){
    rdft((VECTOR_SIZE/WAIT_UNTIL_PERFORM),1,in+i,bitshuffle,weighting);
  }
}

void refresh_buffer_with_incoming_signal(t_sample *buffer, t_sample*in, int n){
  int i=0;
  for(i=0;i<n;i++) {
    buffer[VECTOR_SIZE-n+i]=*in++;
  }
}

t_int *myfft_tilde_perform(t_int *w){

  t_myfft_tilde *x=(t_myfft_tilde *)(w[1]);
  t_sample *in =(t_sample *)(w[2]);
  t_sample *out =(t_sample *)(w[3]);  
  int n =(int)(w[4]);
  

  refresh_buffer_with_incoming_signal(x->circularbuffer,in,n);
  x->bitnumber+= VECTOR_SIZE;
  
  
  if(x->bitnumber>=(WAIT_UNTIL_PERFORM*VECTOR_SIZE)){
    x->bitnumber=0;
    x->complete=1;
    
    apply_blackman(x->circularbuffer,x->window);
    apply_buffered_fft(x->circularbuffer,x->bitshuffle,x->weighting);
  }
  if(x->complete==1){
    copy_to_out_signal(x->circularbuffer,out,n);
  }
  shift_signal(x->circularbuffer,n);
  
  return (w+DSP_ADD_LENGTH+1);
}


void myfft_tilde_dsp(t_myfft_tilde *x, t_signal **sp){
  dsp_add(myfft_tilde_perform, DSP_ADD_LENGTH, x,
    sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}


void myfft_tilde_free(t_myfft_tilde *x){
  free(x->window);
  free(x->bitshuffle);
  free(x->weighting);
  free(x->circularbuffer);
  outlet_free(x->x_out);
}


void *myfft_tilde_new(void)
{
  t_myfft_tilde *m = (t_myfft_tilde *)pd_new(myfft_tilde_class);
  m->x_out = outlet_new(&m->x_obj, &s_signal);
  
  m->window = malloc(VECTOR_SIZE*(sizeof*m->window));
  m->bitshuffle = malloc(2*VECTOR_SIZE*(sizeof*m->bitshuffle));
  m->weighting = malloc(2*VECTOR_SIZE*(sizeof*m->weighting));
  m->circularbuffer = malloc(VECTOR_SIZE*(sizeof*m->circularbuffer));

  m->complete=0;
  m->bitnumber=0;

  init_with_blackman(m->window);
  init_rdft(VECTOR_SIZE, m->bitshuffle, m->weighting);
  
  return (void *)m;
}


void myfft_tilde_setup(void){

  myfft_tilde_class = class_new(gensym("myfft~"),(t_newmethod)myfft_tilde_new,(t_method) myfft_tilde_free,sizeof(t_myfft_tilde),CLASS_DEFAULT,0);
  class_addmethod(myfft_tilde_class,(t_method)myfft_tilde_dsp,gensym("dsp"),0);
  CLASS_MAINSIGNALIN(myfft_tilde_class,t_myfft_tilde,f);
}
