classdef ChangeDimLayer < nnet.layer.Layer
    properties
        Rows  % To store the number of rows
        Cols  % To store the number of columns
    end
    
    methods
        function layer = ChangeDimLayer(name, rows, cols)
            % Layer constructor
            layer.Name = name;
            layer.Type = 'Feature Input';
            layer.Rows = rows;  % Store the number of rows
            layer.Cols = cols;  % Store the number of columns
        end
        
        function Z = predict(layer, X)
            % Flatten the input matrix X into a column vector of size [rows*cols, 1]
            % Use the stored rows and cols to reshape the input X
            Z = reshape(X, layer.Rows, layer.Cols);  % Reshape into a column vector
        end
    end
end
