
library(dplyr)
library(keras)
library(tfdatasets)
library(ggplot2)
boston_housing <- dataset_boston_housing()

c(train_data, train_labels) %<-% boston_housing$train
c(test_data, test_labels) %<-% boston_housing$test



column_names <- c('CRIM', 'ZN', 'INDUS', 'CHAS', 'NOX', 'RM', 'AGE', 
                  'DIS', 'RAD', 'TAX', 'PTRATIO', 'B', 'LSTAT')

train_df <- train_data %>% 
    as_tibble(.name_repair = "minimal") %>% 
    setNames(column_names) %>% 
    mutate(label = train_labels)

test_df <- test_data %>% 
    as_tibble(.name_repair = "minimal") %>% 
    setNames(column_names) %>% 
    mutate(label = test_labels)

spec <- feature_spec(train_df, label ~ . ) %>% 
    step_numeric_column(all_numeric(), normalizer_fn = scaler_standard()) %>% 
    fit()

spec

layer <- layer_dense_features(
    feature_columns = dense_features(spec), 
    dtype = tf$float32
)
layer(train_df)

input <- layer_input_from_dataset(train_df %>% select(-label))

output <- input %>% 
    layer_dense_features(dense_features(spec)) %>% 
    layer_dense(units = 64, activation = "relu") %>%
    layer_dense(units = 64, activation = "relu") %>%
    layer_dense(units = 1) 

model <- keras_model(input, output)

summary(model)

model %>% 
    compile(
        loss = "mse",
        optimizer = optimizer_rmsprop(),
        metrics = list("mean_absolute_error")
    )

build_model <- function() {
    input <- layer_input_from_dataset(train_df %>% select(-label))
    
    output <- input %>% 
        layer_dense_features(dense_features(spec)) %>% 
        layer_dense(units = 64, activation = "relu") %>%
        layer_dense(units = 64, activation = "relu") %>%
        layer_dense(units = 1) 
    
    model <- keras_model(input, output)
    
    model %>% 
        compile(
            loss = "mse",
            optimizer = optimizer_rmsprop(),
            metrics = list("mean_absolute_error")
        )
    
    model
}

print_dot_callback <- callback_lambda(
    on_epoch_end = function(epoch, logs) {
        if (epoch %% 80 == 0) cat("\n")
        cat(".")
    }
)    

model <- build_model()

history <- model %>% fit(
    x = train_df %>% select(-label),
    y = train_df$label,
    epochs = 500,
    validation_split = 0.2,
    verbose = 0,
    callbacks = list(print_dot_callback)
)

plot(history)


model <- keras_model_sequential()
model %>%
    layer_dense(units = c(45, col(TRA.AAf[[1]]))) %>%
    layer_dense(units = 15, activation = "tanh", input_shape = ncol(TRA.AAf[[1]])) %>%
    layer_dense(units = 10, activation = "tanh") %>%
    layer_dense(units = 15, activation = "tanh") %>%
    layer_dense(units = c(45, ncol(TRA.AAf[[1]])))

model %>% compile(
    loss = "mean_squared_error", 
    optimizer = "adam"
)

checkpoint <- callback_model_checkpoint(
    filepath = "model.hdf5", 
    save_best_only = TRUE, 
    period = 1,
    verbose = 1
)

early_stopping <- callback_early_stopping(patience = 5)

x <- tf$Tensor(x.train[[1]], shape = c(45, 28), dtype="int32")
x <- tf$Tensor(x.train, shape=(100, 45, 28), dtype=int32)
x.train <- TRA.AAf[1:100]
x.train <- tf$stack(x.train, axis=0L)


x.test <- TRA.AAf[101:200]
x.test <- tf$stack(x.test, axis=0L)

x_train <- array_reshape(x.train, c(nrow(x.train), 100), order = "F")
x_test <- array_reshape(x.test, c(nrow(x.test), 100), order = "F")


model %>% fit(
    x = x.train, 
    y = x.train, 
    epochs = 100, 
    validation_steps = 100,
    steps_per_epoch = 100,
    batch_size = 32,
    validation_data = list(x.test, x.test), 
    callbacks = list(checkpoint, early_stopping)
)

batch_size <- 100L
original_dim <- 100000L
latent_dim <- 2L
intermediate_dim <- 256L
epochs <- 50L
epsilon_std <- 1.0


x <- layer_input(shape = c(original_dim))
h <- layer_dense(x, intermediate_dim, activation = "relu")
z_mean <- layer_dense(h, latent_dim)
z_log_var <- layer_dense(h, latent_dim)

sampling <- function(arg){
    z_mean <- arg[, 1:(latent_dim)]
    z_log_var <- arg[, (latent_dim + 1):(2 * latent_dim)]
    
    epsilon <- k_random_normal(
        shape = c(k_shape(z_mean)[[1]]), 
        mean=0.,
        stddev=epsilon_std
    )
    
    z_mean + k_exp(z_log_var/2)*epsilon
}

# note that "output_shape" isn't necessary with the TensorFlow backend
z <- layer_concatenate(list(z_mean, z_log_var)) %>% 
    layer_lambda(sampling)

# we instantiate these layers separately so as to reuse them later
decoder_h <- layer_dense(units = intermediate_dim, activation = "relu")
decoder_mean <- layer_dense(units = original_dim, activation = "sigmoid")
h_decoded <- decoder_h(z)
x_decoded_mean <- decoder_mean(h_decoded)

# end-to-end autoencoder
vae <- keras_model(x, x_decoded_mean)

# encoder, from inputs to latent space
encoder <- keras_model(x, z_mean)

# generator, from latent space to reconstructed inputs
decoder_input <- layer_input(shape = latent_dim)
h_decoded_2 <- decoder_h(decoder_input)
x_decoded_mean_2 <- decoder_mean(h_decoded_2)
generator <- keras_model(decoder_input, x_decoded_mean_2)


vae_loss <- function(x, x_decoded_mean){
    xent_loss <- (original_dim/1.0)*loss_binary_crossentropy(x, x_decoded_mean)
    kl_loss <- -0.5*k_mean(1 + z_log_var - k_square(z_mean) - k_exp(z_log_var), axis = -1L)
    xent_loss + kl_loss
}

vae %>% compile(optimizer = "rmsprop", loss = vae_loss)


# Data preparation --------------------------------------------------------

mnist <- keras::dataset_mnist()
x_train <- mnist$train$x/255
x_test <- mnist$test$x/255
x_train <- array_reshape(x_train, c(nrow(x_train), 784), order = "F")
x_test <- array_reshape(x_test, c(nrow(x_test), 784), order = "F")


# Model training ----------------------------------------------------------

vae %>% fit(
    x_train, x_train, 
    shuffle = TRUE, 
    epochs = epochs, 
    batch_size = batch_size, 
    validation_data = list(x_test, x_test)
)
