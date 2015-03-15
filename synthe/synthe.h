#ifndef _synthe_h_
# define _synthe_h_

# include "m_pd.h"

static t_class *synthe_class;

typedef struct _synthe {
  t_object x_obj;

  t_sample f;

  t_inlet *modulatrice;
  t_inlet *messages;
  t_outlet *x_out;

  int *bitshuffle;
  float *weighting;
  float *window;

  t_int bypass;
  t_int autonorm;
  t_float shapeWidth;
} t_synthe;

t_int *synthe_perform(t_int *w);
void synthe_dsp(t_synthe *x, t_signal **sp);
void synthe_free(t_synthe *x);
void *synthe_new(int argc, t_atom *argv);
void synthe_setup(void);
void synthe_messages(t_synthe *x, t_floatarg n, t_floatarg p);

#endif
