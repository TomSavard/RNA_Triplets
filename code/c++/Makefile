# Nom du programme
TARGET = main.exe
TARGET_NUSSINOV = nussinov.exe
TARGET_UTILITIES = utilities.exe



# Compilateur
CC = g++
# Options de compilation
CFLAGS = -Wall -Wextra -std=c++17 -I/opt/anaconda3/envs/INF551_env/include

LDFLAGS = -L/opt/anaconda3/envs/INF551_env/lib -lRNA




# Fichiers source
SRCS = main.cpp global_variables.cpp # Ajoutez ici tous les fichiers source de votre projet
SRCS_NUSSINOV = nussinov.cpp global_variables.cpp utilities.cpp
SRCS_UTILITIES = utilities.cpp global_variables.cpp




# Règle de construction du programme
$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) $(SRCS) -o $(TARGET)

nussinov: $(SRCS_NUSSINOV)
	$(CC) $(CFLAGS) -DNUSSINOV_TEST $(SRCS_NUSSINOV) -o $(TARGET_NUSSINOV) $(LDFLAGS)

utilities: $(SRCS_UTILITIES)
	$(CC) $(CFLAGS) -DUTILITIES_TEST $(SRCS_UTILITIES) -o $(TARGET_UTILITIES)

# Règle pour la construction de tous les fichiers
all: $(TARGET) $(TARGET_NUSSINOV) $(TARGET_UTILITIES)




# Règle pour nettoyer les fichiers objets et l'exécutable
clean:
	rm -f $(TARGET) $(TARGET_NUSSINOV) $(TARGET_UTILITIES)




