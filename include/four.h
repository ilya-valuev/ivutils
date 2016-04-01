# ifdef __cplusplus
extern "C" {
# endif     

typedef double fourtype;
       
void four(fourtype *data, long nn, int isign);
void fold_it(fourtype *data, long nn);
fourtype *four_direct(fourtype *data, long nn, int isign, fourtype w0, fourtype w1,fourtype dw);

# ifdef __cplusplus
}
# endif
