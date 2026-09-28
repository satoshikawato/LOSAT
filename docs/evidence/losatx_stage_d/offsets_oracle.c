/* Comparison only. NCBI blast_gapalign.c3248-3330 owns offsets/failure;
 * blast_stat.c1591-1629 initializes the compiled BLOSUM62 matrix. */
#include <stdio.h>
#include <string.h>
#include <algo/blast/core/blast_gapalign.h>
#include <algo/blast/core/blast_stat.h>
#include <algo/blast/core/blast_encoding.h>
int main(void){
 BlastScoreBlk *sbp=BlastScoreBlkNew(BLASTAA_SEQ_CODE,1);if(!sbp)return 2;
 sbp->name=strdup("BLOSUM62");if(Blast_ScoreBlkMatrixFill(sbp,NULL))return 3;
 Uint1 q[100],s[100];BlastHSP h;Int4 qr,sr;
 for(int mode=0;mode<4;mode++)for(int n=1;n<=51;n++){
 memset(q,1,sizeof(q));memset(s,mode==0?1:21,sizeof(s));memset(&h,0,sizeof(h));h.query.offset=3;h.query.end=3+n;h.subject.offset=7;h.subject.end=7+n+(mode==3?12:0);
 if(mode==2)memset(s+7+n/2,1,n-n/2);
 if(mode==3&&n>11)memset(s+h.subject.end-11,1,11);
 qr=sr=-999;int result=BlastGetOffsetsForGappedAlignment(q,s,sbp,&h,&qr,&sr);
 printf("%d\t%d\t%d\t%d\t%d\n",mode,n,result,qr,sr);
 }
 BlastScoreBlkFree(sbp);return 0;
}
