	  �=  �   k820309    j          14.0        �4_                                                                                                           
       src/runawayelectrongeneration.f RUNAWAYELECTRONGENERATION                      @                             
                            @                              
                        �   @                              
                        �   @                              
                        �   @                             
       %         @                                                  
       #DREICER_GROWTHRATE_CLASSIC%SUM    #DREICER_GROWTHRATE_CLASSIC%SQRT    #DREICER_GROWTHRATE_CLASSIC%ASIN 	   #DREICER_GROWTHRATE_CLASSIC%EXP 
   #Z    #Z0    #NJ    #NSPECIES    #T    #EPAR                                                     SUM                                                  SQRT                                             	     ASIN                                             
     EXP          
                                                          p          5 O p            5 O p                                   
                                                          p          5 O p            5 O p                                   
                                                     
    p          5 O p            5 O p                                    
                                                       
                                      
                
                                      
      %         @                                                 
       #DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%W5    #DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%SUM    #DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%SQRT    #DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%EXP    #DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%RESHAPE    #DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%LOG    #DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%DOT_PRODUCT    #DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%TANH    #DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%MATMUL    #Z    #Z0    #NJ    #NSPECIES    #T    #EPAR                                                                    
                                                     T
W
p          n
      
                 _&�����?        0.3377520000000  n
   
                 C=}���?        0.3740840000000  n
      
                  ~7ݲC��          n
             
                  ��p!׿          n
             
                  M��.�ÿ          n
          
                 �K�Ƽ�?        0.4646470000000  n
   
                 �&2s�K�?        0.6654670000000  n
      
                  �o��ѿ          n
             
                  g�ܶ�ɿ          n
             
                  BZc�	�˿          n
          
                 ��Rb��?        0.2409480000000  n
      
                  :�V�Sҿ          n
             
                  ��>�,Կ          n
             
                  �M~�N�          n
             
                  s�9>Z���          n
             
                  �=���          n
             
                  �:�ҿ          n
             
                  ݕ]0�f�          n
          
                 �6�^��?        0.1843680000000  n
   
                 ���T��?        0.5288490000000  h  p          p          p            p                                                                                                                                         SUM                                                  SQRT                                                  EXP                                                  RESHAPE                                                  LOG                                                  DOT_PRODUCT                                                  TANH                                                  MATMUL          
                                                          p          5 O p            5 O p                                   
                                                          p          5 O p            5 O p                                   
                                                     
    p          5 O p            5 O p                                    
                                                       
                                      
                
                                       
      %         @                               !                   
       #AVALANCHE_GROWTHRATE_CLASSIC%SUM "   #AVALANCHE_GROWTHRATE_CLASSIC%SQRT #   #Z $   #Z0 %   #NJ &   #NSPECIES '   #T (   #EPAR )   #EPS *                                               "     SUM                                             #     SQRT          
                                  $                        p          5 O p            5 O p                                   
                                  %                        p          5 O p            5 O p                                   
                                 &                    
    p          5 O p            5 O p                                    
                                  '                     
                                 (     
                
                                 )     
                
                                 *     
      %         @                               +                   
       #AVALANCHE_GROWTHRATE%SUM ,   #AVALANCHE_GROWTHRATE%SQRT -   #Z .   #Z0 /   #NJ 0   #NSPECIES 1   #T 2   #EPAR 3   #B 4                                               ,     SUM                                             -     SQRT          
                                  .                        p          5 O p            5 O p                                   
                                  /                        p          5 O p            5 O p                                   
                                 0                    
    p          5 O p            5 O p                                    
                                  1                     
                                 2     
                
                                 3     
                
                                 4     
      %         @                                5                   
       #E_CEFF_OVER_E_C%SUM 6   #E_CEFF_OVER_E_C%LOG 7   #E_CEFF_OVER_E_C%SQRT 8   #E_CEFF_OVER_E_C%ABS 9   #Z :   #Z0 ;   #NJ <   #NSPECIES =   #T >   #B ?                                               6     SUM                                             7     LOG                                             8     SQRT                                             9     ABS          
                                  :                        p          5 O p            5 O p                                   
                                  ;                        p          5 O p            5 O p                                   
                                 <                    
    p          5 O p            5 O p                                    
                                  =                     
                                 >     
                
                                 ?     
      %         @                                @                    
       #NE A   #T B             
                                 A     
                
                                 B     
      %         @                                C                    
       #NE D   #T E             
                                 D     
                
                                 E     
      %         @                                F                   
       #P_STAR%SUM G   #P_STAR%SQRT H   #P_STAR%MAX I   #P_STAR%ABS J   #Z K   #Z0 L   #NJ M   #NSPECIES N   #T O   #EPAR P   #B Q                                               G     SUM                                             H     SQRT                                             I     MAX                                             J     ABS          
                                  K                        p          5 O p            5 O p                                   
                                  L                        p          5 O p            5 O p                                   
                                 M                    
 	   p          5 O p            5 O p                                    
                                  N                     
                                 O     
                
                                 P     
                
                                 Q     
      %         @                                R                   
       #NUBAR_D%SUM S   #NUBAR_D%LOG T   #Z U   #Z0 V   #NJ W   #NSPECIES X   #T Y   #P Z                                               S     SUM                                             T     LOG          
                                  U                        p          5 O p            5 O p                                   
                                  V                        p          5 O p            5 O p                                   
                                 W                    
    p          5 O p            5 O p                                    
                                  X                     
                                 Y     
                
                                 Z     
      #         @                                   [                   #CALC_NUBAR_D01%SUM \   #CALC_NUBAR_D01%LOG ]   #Z ^   #Z0 _   #NJ `   #NSPECIES a   #T b   #NUBAR_D0 c   #NUBAR_D1 d                                               \     SUM                                             ]     LOG          
                                  ^                        p          5 O p            5 O p                                   
                                  _                        p          5 O p            5 O p                                   
                                 `                    
    p          5 O p            5 O p                                    
                                  a                     
                                 b     
                                                c     
                                                 d     
       %         @                                e                   
       #NUBAR_S%SUM f   #NUBAR_S%LOG g   #NUBAR_S%SQRT h   #Z i   #Z0 j   #NJ k   #NSPECIES l   #T m   #P n                                               f     SUM                                             g     LOG                                             h     SQRT          
                                  i                        p          5 O p            5 O p                                   
                                  j                        p          5 O p            5 O p                                   
                                 k                    
 	   p          5 O p            5 O p                                    
                                  l                     
                                 m     
                
                                 n     
      #         @                                   o                   #CALC_NUBAR_S01%SUM p   #CALC_NUBAR_S01%LOG q   #Z r   #Z0 s   #NJ t   #NSPECIES u   #T v   #NUBAR_S0 w   #NUBAR_S1 x                                               p     SUM                                             q     LOG          
                                  r                     
   p          5 O p            5 O p                                   
                                  s                        p          5 O p            5 O p                                   
                                 t                    
    p          5 O p            5 O p                                    
                                  u                     
                                 v     
                                                w     
                                                 x     
       %         @                                y                    
       #NE z   #T {             
                                 z     
                
                                 {     
      %         @                                |                   
       #LN_LAMBDA_0%LOG }   #NE ~   #T                                                }     LOG           
                                 ~     
                
                                      
      %         @                                �                   
       #LN_LAMBDA_C%LOG �   #NE �   #T �                                               �     LOG           
                                 �     
                
                                 �     
      %         @                                �                   
       #LN_LAMBDA_EE%LOG �   #LN_LAMBDA_EE%SQRT �   #NE �   #T �   #P �                                               �     LOG                                             �     SQRT           
                                 �     
                
                                 �     
                
                                 �     
      %         @                                �                   
       #LN_LAMBDA_EI%LOG �   #LN_LAMBDA_EI%SQRT �   #NE �   #T �   #P �                                               �     LOG                                             �     SQRT           
                                 �     
                
                                 �     
                
                                 �     
         �   B      fn#fn -   �   @   J   CALCULATE_DREICER_GROWTHRATE /   "  @   J   CALCULATE_AVALANCHE_GROWTHRATE *   b  @   J   CALCULATE_ELECTRIC_FIELDS 0   �  @   J   CALCULATE_COLLISION_FREQUENCIES -   �  @   J   CALCULATE_COULOMB_LOGARITHMS H   "        DREICER_GROWTHRATE_CLASSIC+CALCULATE_DREICER_GROWTHRATE L   :  <      DREICER_GROWTHRATE_CLASSIC%SUM+CALCULATE_DREICER_GROWTHRATE M   v  =      DREICER_GROWTHRATE_CLASSIC%SQRT+CALCULATE_DREICER_GROWTHRATE M   �  =      DREICER_GROWTHRATE_CLASSIC%ASIN+CALCULATE_DREICER_GROWTHRATE L   �  <      DREICER_GROWTHRATE_CLASSIC%EXP+CALCULATE_DREICER_GROWTHRATE J   ,  �   a   DREICER_GROWTHRATE_CLASSIC%Z+CALCULATE_DREICER_GROWTHRATE K   �  �   a   DREICER_GROWTHRATE_CLASSIC%Z0+CALCULATE_DREICER_GROWTHRATE K   t  �   a   DREICER_GROWTHRATE_CLASSIC%NJ+CALCULATE_DREICER_GROWTHRATE Q     @   a   DREICER_GROWTHRATE_CLASSIC%NSPECIES+CALCULATE_DREICER_GROWTHRATE J   X  @   a   DREICER_GROWTHRATE_CLASSIC%T+CALCULATE_DREICER_GROWTHRATE M   �  @   a   DREICER_GROWTHRATE_CLASSIC%EPAR+CALCULATE_DREICER_GROWTHRATE T   �  F      DREICER_GROWTHRATE_CODE_NEURAL_NETWORK+CALCULATE_DREICER_GROWTHRATE W   	  �     DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%W5+CALCULATE_DREICER_GROWTHRATE X   �  <      DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%SUM+CALCULATE_DREICER_GROWTHRATE Y   '  =      DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%SQRT+CALCULATE_DREICER_GROWTHRATE X   d  <      DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%EXP+CALCULATE_DREICER_GROWTHRATE \   �  @      DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%RESHAPE+CALCULATE_DREICER_GROWTHRATE X   �  <      DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%LOG+CALCULATE_DREICER_GROWTHRATE `     D      DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%DOT_PRODUCT+CALCULATE_DREICER_GROWTHRATE Y   `  =      DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%TANH+CALCULATE_DREICER_GROWTHRATE [   �  ?      DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%MATMUL+CALCULATE_DREICER_GROWTHRATE V   �  �   a   DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%Z+CALCULATE_DREICER_GROWTHRATE W   �  �   a   DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%Z0+CALCULATE_DREICER_GROWTHRATE W   $  �   a   DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%NJ+CALCULATE_DREICER_GROWTHRATE ]   �  @   a   DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%NSPECIES+CALCULATE_DREICER_GROWTHRATE V     @   a   DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%T+CALCULATE_DREICER_GROWTHRATE Y   H  @   a   DREICER_GROWTHRATE_CODE_NEURAL_NETWORK%EPAR+CALCULATE_DREICER_GROWTHRATE L   �  �       AVALANCHE_GROWTHRATE_CLASSIC+CALCULATE_AVALANCHE_GROWTHRATE P   d  <      AVALANCHE_GROWTHRATE_CLASSIC%SUM+CALCULATE_AVALANCHE_GROWTHRATE Q   �  =      AVALANCHE_GROWTHRATE_CLASSIC%SQRT+CALCULATE_AVALANCHE_GROWTHRATE N   �  �   a   AVALANCHE_GROWTHRATE_CLASSIC%Z+CALCULATE_AVALANCHE_GROWTHRATE O   �  �   a   AVALANCHE_GROWTHRATE_CLASSIC%Z0+CALCULATE_AVALANCHE_GROWTHRATE O   %  �   a   AVALANCHE_GROWTHRATE_CLASSIC%NJ+CALCULATE_AVALANCHE_GROWTHRATE U   �  @   a   AVALANCHE_GROWTHRATE_CLASSIC%NSPECIES+CALCULATE_AVALANCHE_GROWTHRATE N   	  @   a   AVALANCHE_GROWTHRATE_CLASSIC%T+CALCULATE_AVALANCHE_GROWTHRATE Q   I  @   a   AVALANCHE_GROWTHRATE_CLASSIC%EPAR+CALCULATE_AVALANCHE_GROWTHRATE P   �  @   a   AVALANCHE_GROWTHRATE_CLASSIC%EPS+CALCULATE_AVALANCHE_GROWTHRATE D   �  �       AVALANCHE_GROWTHRATE+CALCULATE_AVALANCHE_GROWTHRATE H   �  <      AVALANCHE_GROWTHRATE%SUM+CALCULATE_AVALANCHE_GROWTHRATE I   �  =      AVALANCHE_GROWTHRATE%SQRT+CALCULATE_AVALANCHE_GROWTHRATE F     �   a   AVALANCHE_GROWTHRATE%Z+CALCULATE_AVALANCHE_GROWTHRATE G   �  �   a   AVALANCHE_GROWTHRATE%Z0+CALCULATE_AVALANCHE_GROWTHRATE G   T  �   a   AVALANCHE_GROWTHRATE%NJ+CALCULATE_AVALANCHE_GROWTHRATE M   �  @   a   AVALANCHE_GROWTHRATE%NSPECIES+CALCULATE_AVALANCHE_GROWTHRATE F   8  @   a   AVALANCHE_GROWTHRATE%T+CALCULATE_AVALANCHE_GROWTHRATE I   x  @   a   AVALANCHE_GROWTHRATE%EPAR+CALCULATE_AVALANCHE_GROWTHRATE F   �  @   a   AVALANCHE_GROWTHRATE%B+CALCULATE_AVALANCHE_GROWTHRATE :   �  �       E_CEFF_OVER_E_C+CALCULATE_ELECTRIC_FIELDS >   �  <      E_CEFF_OVER_E_C%SUM+CALCULATE_ELECTRIC_FIELDS >     <      E_CEFF_OVER_E_C%LOG+CALCULATE_ELECTRIC_FIELDS ?   X  =      E_CEFF_OVER_E_C%SQRT+CALCULATE_ELECTRIC_FIELDS >   �  <      E_CEFF_OVER_E_C%ABS+CALCULATE_ELECTRIC_FIELDS <   �  �   a   E_CEFF_OVER_E_C%Z+CALCULATE_ELECTRIC_FIELDS =   u  �   a   E_CEFF_OVER_E_C%Z0+CALCULATE_ELECTRIC_FIELDS =     �   a   E_CEFF_OVER_E_C%NJ+CALCULATE_ELECTRIC_FIELDS C   �  @   a   E_CEFF_OVER_E_C%NSPECIES+CALCULATE_ELECTRIC_FIELDS <   �  @   a   E_CEFF_OVER_E_C%T+CALCULATE_ELECTRIC_FIELDS <   =   @   a   E_CEFF_OVER_E_C%B+CALCULATE_ELECTRIC_FIELDS .   }   _       E_C+CALCULATE_ELECTRIC_FIELDS 1   �   @   a   E_C%NE+CALCULATE_ELECTRIC_FIELDS 0   !  @   a   E_C%T+CALCULATE_ELECTRIC_FIELDS .   \!  _       E_D+CALCULATE_ELECTRIC_FIELDS 1   �!  @   a   E_D%NE+CALCULATE_ELECTRIC_FIELDS 0   �!  @   a   E_D%T+CALCULATE_ELECTRIC_FIELDS 6   ;"  �       P_STAR+CALCULATE_AVALANCHE_GROWTHRATE :   	#  <      P_STAR%SUM+CALCULATE_AVALANCHE_GROWTHRATE ;   E#  =      P_STAR%SQRT+CALCULATE_AVALANCHE_GROWTHRATE :   �#  <      P_STAR%MAX+CALCULATE_AVALANCHE_GROWTHRATE :   �#  <      P_STAR%ABS+CALCULATE_AVALANCHE_GROWTHRATE 8   �#  �   a   P_STAR%Z+CALCULATE_AVALANCHE_GROWTHRATE 9   �$  �   a   P_STAR%Z0+CALCULATE_AVALANCHE_GROWTHRATE 9   B%  �   a   P_STAR%NJ+CALCULATE_AVALANCHE_GROWTHRATE ?   �%  @   a   P_STAR%NSPECIES+CALCULATE_AVALANCHE_GROWTHRATE 8   &&  @   a   P_STAR%T+CALCULATE_AVALANCHE_GROWTHRATE ;   f&  @   a   P_STAR%EPAR+CALCULATE_AVALANCHE_GROWTHRATE 8   �&  @   a   P_STAR%B+CALCULATE_AVALANCHE_GROWTHRATE 8   �&  �       NUBAR_D+CALCULATE_COLLISION_FREQUENCIES <   �'  <      NUBAR_D%SUM+CALCULATE_COLLISION_FREQUENCIES <   �'  <      NUBAR_D%LOG+CALCULATE_COLLISION_FREQUENCIES :   (  �   a   NUBAR_D%Z+CALCULATE_COLLISION_FREQUENCIES ;   �(  �   a   NUBAR_D%Z0+CALCULATE_COLLISION_FREQUENCIES ;   K)  �   a   NUBAR_D%NJ+CALCULATE_COLLISION_FREQUENCIES A   �)  @   a   NUBAR_D%NSPECIES+CALCULATE_COLLISION_FREQUENCIES :   /*  @   a   NUBAR_D%T+CALCULATE_COLLISION_FREQUENCIES :   o*  @   a   NUBAR_D%P+CALCULATE_COLLISION_FREQUENCIES ?   �*  �       CALC_NUBAR_D01+CALCULATE_COLLISION_FREQUENCIES C   o+  <      CALC_NUBAR_D01%SUM+CALCULATE_COLLISION_FREQUENCIES C   �+  <      CALC_NUBAR_D01%LOG+CALCULATE_COLLISION_FREQUENCIES A   �+  �   a   CALC_NUBAR_D01%Z+CALCULATE_COLLISION_FREQUENCIES B   �,  �   a   CALC_NUBAR_D01%Z0+CALCULATE_COLLISION_FREQUENCIES B   /-  �   a   CALC_NUBAR_D01%NJ+CALCULATE_COLLISION_FREQUENCIES H   �-  @   a   CALC_NUBAR_D01%NSPECIES+CALCULATE_COLLISION_FREQUENCIES A   .  @   a   CALC_NUBAR_D01%T+CALCULATE_COLLISION_FREQUENCIES H   S.  @   a   CALC_NUBAR_D01%NUBAR_D0+CALCULATE_COLLISION_FREQUENCIES H   �.  @   a   CALC_NUBAR_D01%NUBAR_D1+CALCULATE_COLLISION_FREQUENCIES 8   �.  �       NUBAR_S+CALCULATE_COLLISION_FREQUENCIES <   �/  <      NUBAR_S%SUM+CALCULATE_COLLISION_FREQUENCIES <   �/  <      NUBAR_S%LOG+CALCULATE_COLLISION_FREQUENCIES =   0  =      NUBAR_S%SQRT+CALCULATE_COLLISION_FREQUENCIES :   ?0  �   a   NUBAR_S%Z+CALCULATE_COLLISION_FREQUENCIES ;   �0  �   a   NUBAR_S%Z0+CALCULATE_COLLISION_FREQUENCIES ;   �1  �   a   NUBAR_S%NJ+CALCULATE_COLLISION_FREQUENCIES A   +2  @   a   NUBAR_S%NSPECIES+CALCULATE_COLLISION_FREQUENCIES :   k2  @   a   NUBAR_S%T+CALCULATE_COLLISION_FREQUENCIES :   �2  @   a   NUBAR_S%P+CALCULATE_COLLISION_FREQUENCIES ?   �2  �       CALC_NUBAR_S01+CALCULATE_COLLISION_FREQUENCIES C   �3  <      CALC_NUBAR_S01%SUM+CALCULATE_COLLISION_FREQUENCIES C   �3  <      CALC_NUBAR_S01%LOG+CALCULATE_COLLISION_FREQUENCIES A   #4  �   a   CALC_NUBAR_S01%Z+CALCULATE_COLLISION_FREQUENCIES B   �4  �   a   CALC_NUBAR_S01%Z0+CALCULATE_COLLISION_FREQUENCIES B   k5  �   a   CALC_NUBAR_S01%NJ+CALCULATE_COLLISION_FREQUENCIES H   6  @   a   CALC_NUBAR_S01%NSPECIES+CALCULATE_COLLISION_FREQUENCIES A   O6  @   a   CALC_NUBAR_S01%T+CALCULATE_COLLISION_FREQUENCIES H   �6  @   a   CALC_NUBAR_S01%NUBAR_S0+CALCULATE_COLLISION_FREQUENCIES H   �6  @   a   CALC_NUBAR_S01%NUBAR_S1+CALCULATE_COLLISION_FREQUENCIES 6   7  _       NU_EE+CALCULATE_COLLISION_FREQUENCIES 9   n7  @   a   NU_EE%NE+CALCULATE_COLLISION_FREQUENCIES 8   �7  @   a   NU_EE%T+CALCULATE_COLLISION_FREQUENCIES 9   �7  t       LN_LAMBDA_0+CALCULATE_COULOMB_LOGARITHMS =   b8  <      LN_LAMBDA_0%LOG+CALCULATE_COULOMB_LOGARITHMS <   �8  @   a   LN_LAMBDA_0%NE+CALCULATE_COULOMB_LOGARITHMS ;   �8  @   a   LN_LAMBDA_0%T+CALCULATE_COULOMB_LOGARITHMS 9   9  t       LN_LAMBDA_C+CALCULATE_COULOMB_LOGARITHMS =   �9  <      LN_LAMBDA_C%LOG+CALCULATE_COULOMB_LOGARITHMS <   �9  @   a   LN_LAMBDA_C%NE+CALCULATE_COULOMB_LOGARITHMS ;   :  @   a   LN_LAMBDA_C%T+CALCULATE_COULOMB_LOGARITHMS :   N:  �       LN_LAMBDA_EE+CALCULATE_COULOMB_LOGARITHMS >   �:  <      LN_LAMBDA_EE%LOG+CALCULATE_COULOMB_LOGARITHMS ?   ;  =      LN_LAMBDA_EE%SQRT+CALCULATE_COULOMB_LOGARITHMS =   Z;  @   a   LN_LAMBDA_EE%NE+CALCULATE_COULOMB_LOGARITHMS <   �;  @   a   LN_LAMBDA_EE%T+CALCULATE_COULOMB_LOGARITHMS <   �;  @   a   LN_LAMBDA_EE%P+CALCULATE_COULOMB_LOGARITHMS :   <  �       LN_LAMBDA_EI+CALCULATE_COULOMB_LOGARITHMS >   �<  <      LN_LAMBDA_EI%LOG+CALCULATE_COULOMB_LOGARITHMS ?   �<  =      LN_LAMBDA_EI%SQRT+CALCULATE_COULOMB_LOGARITHMS =   &=  @   a   LN_LAMBDA_EI%NE+CALCULATE_COULOMB_LOGARITHMS <   f=  @   a   LN_LAMBDA_EI%T+CALCULATE_COULOMB_LOGARITHMS <   �=  @   a   LN_LAMBDA_EI%P+CALCULATE_COULOMB_LOGARITHMS 