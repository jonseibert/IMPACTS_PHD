function Lv = compute_Lv(T,Lv_chart,T_chart)

Lv = zeros(size(T));
chart_idx_a = 0;
chart_idx_b = 0;
for z_idx = 1:length(T)
    T_cur = T(z_idx);
    for lv_idx = 1:length(Lv_chart)
        T_chart_cur = T_chart(lv_idx);
        chart_idx_a = lv_idx;
        if(T_chart_cur > T_cur)
            if(lv_idx == 1)
                T_chart_prev = T_chart_cur;
                chart_idx_b = lv_idx;
            else
                T_chart_prev = T_chart(lv_idx-1);
                chart_idx_b = lv_idx-1;
            end
            break;
        elseif(lv_idx == length(Lv_chart))
            T_chart_prev = T_chart_cur;
            chart_idx_b = lv_idx;
        end
    end
    Lv(z_idx) = linear_interpolate(T_chart_cur,T_chart_prev,T_cur,Lv_chart(chart_idx_a),Lv_chart(chart_idx_b));
end