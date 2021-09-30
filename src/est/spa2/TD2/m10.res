Thu Sep 30 07:53:11 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	Two-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa2/TD2/dat10.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
$PK
ET1 = EXP(ETA(1)*THETA(6))
ET2 = EXP(ETA(2)*THETA(7))
ET3 = EXP(ETA(3)*THETA(8))
ET4 = EXP(ETA(4)*THETA(9))
ET5 = EXP(ETA(5)*THETA(10))

CL = 5.0 * THETA(1) * ET1
V2 = 35  * THETA(2) * ET2
Q  = 50  * THETA(3) * ET3
V3 = 50  * THETA(4) * ET4
KA = 0.7 * THETA(5) * ET5
SC = V2
$ERROR
CVERR = 0.05
W = THETA(11)*F*CVERR

Y 	= F + W*ERR(1)

$THETA
(0,1) ; CL
(0,1) ; V2
(0,1) ; Q
(0,1) ; V3
(0,1) ; KA
(0,1) ; IIVCL
(0,1) ; IIVV2
(0,1) ; IIVQ
(0,1) ; IIVV3
(0,1) ; IIVKA
(0,1) ; CVPropErr

$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       30 SEP 2021
Days until program expires : 199
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 template control stream
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      700
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      600
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:  11
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 TWO COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN4)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K23)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K32)
   BASIC PK PARAMETER NO.  5: ABSORPTION RATE (KA)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V2, Q, V3 TO K, K23, K32 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         PERIPH.      ON         NO         YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            6           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
1


 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            10000
 NO. OF SIG. FIGURES REQUIRED:            2
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): m10.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2415.54607521987        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6415E+02  3.6701E+01 -2.9164E+01  6.6768E+01  1.1986E+01  2.5074E+01  2.0388E+01  2.0029E+00 -1.0736E+01  2.1343E+01
            -7.0146E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2423.76537797764        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.8222E-01  1.0301E+00  1.1282E+00  1.0000E+00  1.0581E+00  1.0991E+00  9.2593E-01  9.9735E-01  1.0553E+00  9.3236E-01
             1.0613E+00
 PARAMETER:  8.2065E-02  1.2967E-01  2.2063E-01  1.0005E-01  1.5649E-01  1.9446E-01  2.3048E-02  9.7349E-02  1.5387E-01  2.9968E-02
             1.5947E-01
 GRADIENT:  -1.6593E+00  2.8638E+00  3.2108E+00  4.1243E+00 -9.3324E+00  8.1105E+00  9.2637E+00 -6.9898E+00 -2.0995E+00 -2.6921E+00
            -1.8056E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2426.47511527258        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  9.8508E-01  1.1213E+00  1.1621E+00  9.5714E-01  1.1414E+00  1.0927E+00  6.7247E-01  1.2626E+00  1.1396E+00  1.0488E+00
             1.0506E+00
 PARAMETER:  8.4968E-02  2.1448E-01  2.5027E-01  5.6195E-02  2.3228E-01  1.8865E-01 -2.9680E-01  3.3316E-01  2.3070E-01  1.4767E-01
             1.4935E-01
 GRADIENT:   3.3012E+00  8.3211E+00 -1.0185E+01  2.4336E+01  6.8260E+00  6.0715E+00  2.1605E+00  9.2278E-01 -2.7033E-01  3.8306E+00
            -8.3370E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2427.10357150848        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      529
 NPARAMETR:  9.8336E-01  1.1400E+00  1.2128E+00  9.2993E-01  1.1645E+00  1.0747E+00  5.9023E-01  1.3314E+00  1.1698E+00  1.0566E+00
             1.0594E+00
 PARAMETER:  8.3215E-02  2.3102E-01  2.9292E-01  2.7356E-02  2.5226E-01  1.7201E-01 -4.2725E-01  3.8620E-01  2.5681E-01  1.5508E-01
             1.5771E-01
 GRADIENT:   4.5585E-02  1.0388E+00  1.5289E-01 -9.8745E-04 -1.8311E+00 -2.6785E-01  6.5604E-01 -2.4613E-01 -1.8659E-02  3.6181E-01
             1.3358E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2427.22708402426        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      711
 NPARAMETR:  9.8373E-01  1.1187E+00  1.3351E+00  9.4334E-01  1.1925E+00  1.0756E+00  4.1858E-01  1.5532E+00  1.1869E+00  1.0680E+00
             1.0584E+00
 PARAMETER:  8.3598E-02  2.1214E-01  3.8904E-01  4.1676E-02  2.7608E-01  1.7289E-01 -7.7089E-01  5.4030E-01  2.7132E-01  1.6581E-01
             1.5678E-01
 GRADIENT:   1.1845E+00  2.6529E+00  9.6960E-01  6.0555E-01 -1.3957E+00  4.0697E-01  3.5461E-02 -8.5150E-01 -1.5573E+00 -8.6854E-01
            -3.4914E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2427.37122567720        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      889
 NPARAMETR:  9.8320E-01  1.1138E+00  1.3912E+00  9.4207E-01  1.2138E+00  1.0727E+00  1.7201E-01  1.7137E+00  1.2317E+00  1.0986E+00
             1.0570E+00
 PARAMETER:  8.3052E-02  2.0774E-01  4.3019E-01  4.0326E-02  2.9375E-01  1.7017E-01 -1.6602E+00  6.3863E-01  3.0838E-01  1.9402E-01
             1.5544E-01
 GRADIENT:   4.0610E-01  7.0807E-01 -1.5699E+00  1.3400E+00  2.5070E+00 -5.6467E-01  1.1740E-01  4.8503E-01 -1.1338E+00 -3.1001E-01
             3.8311E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2427.43823426748        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1065
 NPARAMETR:  9.8259E-01  1.1119E+00  1.4097E+00  9.4087E-01  1.2155E+00  1.0747E+00  4.5592E-02  1.7360E+00  1.2459E+00  1.1068E+00
             1.0574E+00
 PARAMETER:  8.2434E-02  2.0611E-01  4.4340E-01  3.9053E-02  2.9520E-01  1.7202E-01 -2.9880E+00  6.5161E-01  3.1988E-01  2.0148E-01
             1.5585E-01
 GRADIENT:  -6.8039E-01  2.2370E-01  4.9770E-01 -3.1046E-01 -6.9254E-01  2.0736E-01  2.1303E-02 -4.9353E-02  3.7034E-01  1.5734E-01
             1.2212E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2427.44590238603        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1240
 NPARAMETR:  9.8295E-01  1.1104E+00  1.4085E+00  9.4192E-01  1.2151E+00  1.0740E+00  1.0440E-02  1.7374E+00  1.2443E+00  1.1061E+00
             1.0564E+00
 PARAMETER:  8.2803E-02  2.0475E-01  4.4252E-01  4.0170E-02  2.9480E-01  1.7141E-01 -4.4621E+00  6.5237E-01  3.1861E-01  2.0084E-01
             1.5484E-01
 GRADIENT:   4.2902E-02 -4.4381E-02 -2.7158E-02 -7.2909E-02  3.6476E-02 -4.0520E-02  2.4621E-03  1.1948E-02  5.1353E-02  2.2593E-02
             1.3991E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2427.44681771460        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1406
 NPARAMETR:  9.8398E-01  1.1104E+00  1.4087E+00  9.4203E-01  1.2150E+00  1.0760E+00  1.0000E-02  1.7378E+00  1.2444E+00  1.1060E+00
             1.0562E+00
 PARAMETER:  8.3850E-02  2.0474E-01  4.4266E-01  4.0286E-02  2.9471E-01  1.7322E-01 -4.6176E+00  6.5262E-01  3.1864E-01  2.0078E-01
             1.5471E-01
 GRADIENT:   2.0497E+00  1.5161E-01  2.8341E-02  6.3875E-02 -1.3793E-01  6.8479E-01  0.0000E+00  2.1599E-02  6.4521E-02  1.5056E-02
             8.1845E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1406
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.5998E-04 -1.0773E-03 -2.8480E-02 -5.6193E-04 -2.0538E-02
 SE:             2.9894E-02  3.1379E-04  1.8401E-02  2.9468E-02  2.5188E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8505E-01  5.9664E-04  1.2169E-01  9.8479E-01  4.1484E-01

 ETASHRINKSD(%)  1.0000E-10  9.8949E+01  3.8353E+01  1.2784E+00  1.5619E+01
 ETASHRINKVR(%)  1.0000E-10  9.9989E+01  6.1997E+01  2.5405E+00  2.8798E+01
 EBVSHRINKSD(%)  3.1118E-01  9.9123E+01  4.0588E+01  1.7441E+00  1.4461E+01
 EBVSHRINKVR(%)  6.2140E-01  9.9992E+01  6.4702E+01  3.4579E+00  2.6831E+01
 RELATIVEINF(%)  9.9357E+01  1.5945E-03  1.8407E+01  2.7114E+01  3.1190E+01
 EPSSHRINKSD(%)  2.9630E+01
 EPSSHRINKVR(%)  5.0480E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2427.4468177146036     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1324.7205778689965     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.72
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.59
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2427.447       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.84E-01  1.11E+00  1.41E+00  9.42E-01  1.21E+00  1.08E+00  1.00E-02  1.74E+00  1.24E+00  1.11E+00  1.06E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  9.00E-02
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  3.00E-01
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.00E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        9.85E+02
 
 TH 2
+       -8.49E+00  7.41E+02
 
 TH 3
+        9.55E-01  8.77E+01  6.85E+01
 
 TH 4
+       -6.23E+00  5.38E+02 -1.13E+01  7.81E+02
 
 TH 5
+       -7.96E-02 -2.81E+02 -1.08E+02 -2.11E+01  4.51E+02
 
 TH 6
+        1.89E-01 -7.72E-01  2.11E-01 -1.33E+00 -1.36E-01  1.70E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        1.63E-01 -1.07E+01 -1.37E+01 -1.92E+00 -4.52E+00 -8.65E-04  0.00E+00  1.61E+01
 
 TH 9
+        1.43E+00 -1.08E+02  1.38E+00  4.11E+00  1.97E+00 -7.39E-01  0.00E+00 -4.86E-01  1.20E+02
 
 TH10
+        4.17E-01 -3.77E+01 -8.11E-01  1.58E+00 -2.06E+01  2.77E-01  0.00E+00  3.83E+00  1.42E+00  9.68E+01
 
 TH11
+       -6.46E+00 -2.62E+01 -3.86E+00 -9.89E+00 -2.68E+00  2.27E+00  0.00E+00  1.21E+01  6.23E+00  1.19E+01  4.88E+02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       34.383
Stop Time:
Thu Sep 30 07:53:47 CDT 2021
