function Meta_Data=mod_som_get_filters_name_efe(Meta_Data)


switch Meta_Data.Firmware.version
    case{'MADRE2.1','SOM'}
        H.shear = 'sinc4';
        H.FPO7  = 'sinc4';
        Meta_Data.epsi.s1.ADCfilter=H.shear;
        Meta_Data.epsi.s2.ADCfilter=H.shear;
        Meta_Data.epsi.t1.ADCfilter=H.FPO7;
        Meta_Data.epsi.t2.ADCfilter=H.FPO7;

end