#include <stdlib.h>
#include <math.h>
#include "myfft_fft.c"
#include "synthe.h"

#define PI 3.14159265358979323846
#define DSP_ADD_LENGTH 5
#define VECTOR_SIZE 2048

/* Fill vector with begin values according to BlackMan
 * @param A vector to fill
 */
void init_with_blackman(t_sample *in){
  int i;
  for (i=0;i<VECTOR_SIZE;i++) {
    in[i] = (float) (0.42-0.5*(cos (2*PI * i/VECTOR_SIZE) )+0.08*(cos (4*PI*i/VECTOR_SIZE)) );
  }
}

/* Copy a signal into another
 * @param A vector model
 * @param A vector copied from the model
 * @param Size of vector 
 */
void copy_to_out_signal(t_sample *in, t_sample *out, int n){
  int i;
  for(i=0;i<n;i++){
    out[i]= in[i];
  }  
}

/* Naive euclidienne distance
 * @param first element for the distance
 * @param second element for the distance
 */
float distance_euclidienne(float x, float y){
  return sqrtf(x*x+y*y);
}

/* Apply BlackMan to a vector
 * @param A vector to be applied
 * @param A vector with blackman values
 */
void apply_blackman(t_sample *in1, t_sample *in2,  t_sample *blackman, int n){
  int i;
  for(i=0;i<n;i++){
    in1[i] = blackman[i] * in2[i];
  }
}

/* Trigger function for Pure Data DSP 
 *  @param  Vector which contains DSP add informations 
 *  @return Size of Vector and add one 
 */
t_int *synthe_perform(t_int *w){

  t_synthe *x=(t_synthe*)(w[1]);  
  t_sample *in1=(t_sample *)(w[2]);  
  t_sample *in2=(t_sample *)(w[3]);  
  t_sample *out=(t_sample *)(w[4]);  
  int n=(int)(w[5]),i=0,j=0;

  t_sample *dup1,*dup2,*tmp1,*tmp2,*res;
  t_sample a1,a2,b1,b2;
  float ampSum,freqSum,factor;
	
  dup1=malloc(n*(sizeof*dup1));
  dup2=malloc(n*(sizeof*dup2));
  tmp1=malloc(n*(sizeof*tmp1));
  tmp2=malloc(n*(sizeof*tmp2));
  res=malloc(n*(sizeof*res));
  
  // Recopier l'in1 dans l'out et sortir
  if (x->bypass==1) {
    copy_to_out_signal(in1,out,n);
    return (w+6);
  }
  
  // Appliquer la fenetre aux buffer dupliqués (dup1 et dup2)
  apply_blackman(dup1,in1,x->window,n);
  apply_blackman(dup2,in2,x->window,n);
  
  // Application de la fonction rfdt aux deux buffers
  rdft(n,1,dup1,x->bitshuffle,x->weighting);
  rdft(n,1,dup2,x->bitshuffle,x->weighting);
	
  // Conversion des valeurs complexe en coordonne polaires
  for(i=0;i<n;i=i+2){
    a1=dup1[i];b1=dup1[i+1];
    a2=dup2[i];b2=dup2[i+1];
    tmp1[i]=distance_euclidienne(a1,b1);
    tmp1[i+1]=-atan2(b1,a1);
    tmp2[i]=distance_euclidienne(a2,b2);
    tmp2[i+1]=-atan2(b2,a2);
  }
  for(i=0;i<n;i=(i+x->shapeWidth*2)) {
    ampSum=0;
    freqSum=0;
    for(j=0;j<(x->shapeWidth)*2;j=j+2) {
      ampSum+=tmp2[i+j];
      freqSum+=tmp1[i+j];
    }
    factor=ampSum/freqSum;
    for(j=0;j<(x->shapeWidth)*2;j=j+2) {
      tmp1[i+j]*=factor;
    }
  }

  // reconvertir en valeurs complexes
  for(i=0;i<n;i=i+2){
    res[i]=tmp1[i]*cos(tmp1[i+1]);
    res[i+1]=-tmp1[i]*sin(tmp1[i+1]);
  }
  rdft(n,-1,res,x->bitshuffle,x->weighting);	
  if(x->autonorm==1){
    for(i=0;i<n;i++){
      res[i]=res[i]/500;
    }
  }
  copy_to_out_signal(res,out,n);
	
  free(dup1);
  free(dup2);
  free(tmp1);
  free(tmp2);
  free(res);
	
  return w+6;
}

/* Initialize synthe on PureData 
   @param (autonorm value)? (bypass value)?
 */
void *synthe_new(int argc, t_atom *argv){

  t_synthe *m = (t_synthe *)pd_new(synthe_class);

  m->window = malloc(VECTOR_SIZE * (sizeof * m->window));
  m->bitshuffle = malloc(VECTOR_SIZE * 2 * (sizeof * m->bitshuffle));
  m->weighting = malloc(VECTOR_SIZE * 2 * (sizeof * m->weighting));
  init_with_blackman(m->window);
  init_rdft(VECTOR_SIZE, m->bitshuffle, m->weighting);

  switch (argc) {
  case 1 :
    m->autonorm = atom_getint(argv);
    m->bypass = 1;
    m->shapeWidth = 0.5;
  case 2 :
    m->bypass = atom_getint(argv + (1 * sizeof(t_atom)));
    m->autonorm = 1;
    m->shapeWidth = 0.5;
  default :
    m->autonorm = 1;
    m->bypass = 1;
    m->shapeWidth = 0.5;
  }
  m->modulatrice = inlet_new(&m->x_obj, &m->x_obj.ob_pd,&s_signal, &s_signal);
  floatinlet_new(&m->x_obj, &m->shapeWidth);
  m->messages = inlet_new(&m->x_obj,&m->x_obj.ob_pd,gensym("list"),gensym("messages"));
  m->x_out = outlet_new(&m->x_obj, &s_signal);

  	
  return (void*)m;
}

/* Add myfft~ into DSP stack 
 *  @param myfft~ object
 *  @param input signal
 */
void synthe_dsp(t_synthe *x, t_signal **sp){
  dsp_add(synthe_perform, DSP_ADD_LENGTH, x,  
          sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[0]->s_n); 
}

/* Message inlet to manage autonorm and bypass 
 * @param autonorm value
 * @param bypass value
*/
void synthe_messages(t_synthe *x, t_floatarg n, t_floatarg p){
  x->autonorm=n;
  x->bypass=p;
}

/* Free the memory of myfft~ 
 *  @param myfft~ object 
 */
void synthe_free(t_synthe *x){
  free(x->window);
  free(x->bitshuffle);
  free(x->weighting);

  inlet_free(x->modulatrice);
  inlet_free(x->messages);
  outlet_free(x->x_out);
}

/* Setup synthe on PureData */
void synthe_setup(void){
  synthe_class = class_new(gensym("synthe"),(t_newmethod)synthe_new,(t_method) synthe_free, sizeof(t_synthe),CLASS_DEFAULT, A_GIMME, 0);
  class_addmethod(synthe_class, (t_method)synthe_dsp, gensym("dsp"), 0);
  class_addmethod(synthe_class, (t_method)synthe_messages,gensym("messages"), A_DEFFLOAT, A_DEFFLOAT, 0);
  CLASS_MAINSIGNALIN(synthe_class, t_synthe, f);
}
