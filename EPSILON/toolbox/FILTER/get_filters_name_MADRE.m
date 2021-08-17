function Meta_Data=get_filters_name_MADRE(Meta_Data)


switch Meta_Data.Firmware.version
    case{'MADRE2.1'}
        H.shear = 'sinc4';
        H.FPO7  = 'sinc4';
        Meta_Data.epsi.s1.ADCfilter=H.shear;
        Meta_Data.epsi.s2.ADCfilter=H.shear;
        Meta_Data.epsi.t1.ADCfilter=H.FPO7;
        Meta_Data.epsi.t2.ADCfilter=H.FPO7;

end