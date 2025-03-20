function err = RMSE(A,B)
    data_dif = A-B;
    data_dif_clean = data_dif(~isnan(data_dif));
    err = sqrt(sum(data_dif_clean.^2)/(numel(data_dif_clean)));
end