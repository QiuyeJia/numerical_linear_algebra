function getrfRook!(A)
    # Rook pivoting for rank-revealing LU
    function getrfRook!(A0)
    # Rook pivoting for rank-revealing LU
    n = size(A0,1)
    P_row = collect(1:n); P_col = collect(1:n);
    for k=1:r
        # Rook pivoting
        row = 1; row0 = 0; col = 1; col0 = 0
        while row != row0 || col != col0
            row0, col0 = row, col # Save old values
            row_A0 = abs.(A0[row+k-1, k:end]) # Search in pivots' row
            col = argmax(row_A0)
            col_A0 = abs.(A0[k:end, col+k-1]) # Search in pivot's column
            row = argmax(col_A0)
        end
        # If we reach this line, this means that the pivot is the largest
        # in its row and column.
        row += k-1; col += k-1
        # Swap rows and columns
        P_row[k], P_row[row] = P_row[row], P_row[k]
        P_col[k], P_col[col] = P_col[col], P_col[k]
        for j=1:n
            A0[k,j],A0[row,j] = A0[row,j],A0[k,j]
        end
        for i=1:n
            A0[i,k],A0[i,col] = A0[i,col],A0[i,k]
        end

        # Perform usual LU step
        if A0[k,k] != 0
            for i=k+1:r
                A0[i,k] /= A0[k,k]
            end
            for j=k+1:r, i=k+1:r
                A0[i,j] -= A0[i,k] * A0[k,j]
            end
        end
    end
    return P_row, P_col
end
end
