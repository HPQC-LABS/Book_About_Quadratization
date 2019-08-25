function pred = predict(x)
    global data
    a1 = relu(x*data.W1 + data.b1);
    a2 = relu(a1*data.W2 + data.b2);
    a3 = relu(a2*data.W3 + data.b3);
    a4 = sigmoid(a3*data.W4 + data.b4);
    
    if a4 >= 0.5
        pred = 1;
    else
        pred = 0;
    end
end
