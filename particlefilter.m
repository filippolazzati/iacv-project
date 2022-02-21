

while nframe < size(frames, 4)
    frame = frames(:,:,:,nframe);
    % ...

    X = F * X;
    X(1:2,:) = X(1:2,:) + Xstd_pos * randn(2, N);
    X(3:4,:) = X(3:4,:) + Xstd_vel * randn(2, N);

    avg_pos_x = mean(X(1,:));
    avg_pos_y = mean(X(2,:));

    selected_cc = nan(1, 2);
    if size(cc, 1) > 0
        for c = 1:size(cc, 1)
            if (cc(c, 1) > avg_pos_x - window_width && cc(c, 1) < avg_pos_x + window_width) || (cc(c, 2) > avg_pos_y - window_height && cc(c, 2) < avg_pos_y + window_height)
                selected_cc = cc(c, :);
                break;
            end
        end
    end

    if ~isnan(selected_cc(1))
        for p = 1:n_particles
            w(p) = abs(selected_cc(1) - X(1, p)) + abs(selected_cc(2) - X(2, p)) + randn * Xstd_measure;
        end
        w = w ./ sum(w);
        pdf = cumsum(w);
        [~, ~, I] = histcounts(rand(1, n_particles), pdf);
        X = X(:, I + 1); % +1 to avoid zero index
    end

    predicted_x = mean(X(1,:));
    predicted_y = mean(X(2,:));

    disp([predicted_x, predicted_y]);
end









