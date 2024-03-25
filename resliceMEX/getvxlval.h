#ifndef  GETVXLVAL_H
#define  GETVXLVAL_H

int index(const int *ps, int i,int j, int k);
void ind2sub(const int *ps, const int ind, int &i, int &j, int &k);
bool inbound(const double *pvxl, const int *ps, int interpm, int pad);
bool getvxlval(const double *pvxlo, const int *ps, const double *pimg, const double *pmv2v, int interpm, int padding, double &val);
double getmin(int isize, const double *pval);
double getmax(int isize, const double *pval);
int binind1d(double dbr, double dmin, int numbins, double val, double thmin, double thmax);
double gethist(int isize, const double *pvxl, double dmin, double dmax, int numbins, double *pval, double *pbins, double thmin, double thmax);
bool binind2d(int numbins1, int numbins2, double dbr1, double dbr2, double dmin1, double dmin2, float val1, float val2, double thmin1, double thmax1, double thmin2, double thmax2, int* pind);
double getjhist(int vsize, const double *pvxlval1, const double *pvxlval2, double dmin1, double dmax1, int numbins1, double thmin1 ,double thmax1, double dmin2, double dmax2, int numbins2, double thmin2 ,double thmax2, double *pval, double *pbins1, double *pbins2, double *pslog10);

#endif
