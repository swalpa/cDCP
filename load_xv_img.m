function d = load_xv_img(fn)
    
    fid = fopen(fn,'r');
    data = fread(fid);
    fclose(fid);

    d = data(1025:end);
    d = reshape(d, 64,64);
end