function I_3=meafilt2(I_noise)
I_3=fspecial('average',[3,3]);%3*3meanfilter
I_3=imfilter(I_noise,I_3);