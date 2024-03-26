#ifndef Pythia8_Win32Defs_H
#define Pythia8_Win32Defs_H

#include <io.h>
#include <Windows.h>
#include <Psapi.h>
#ifndef RTLD_DEFAULT
#define RTLD_DEFAULT ((void *)::GetModuleHandle(NULL))
#define dlopen(library_name, flags) ::LoadLibrary(library_name)
#define dlclose(library) ::FreeLibrary((HMODULE)library)
#endif // RTLD_DEFAULT

char *dlerror();
FARPROC dlsym(void *, const char *);

#endif // Pythia8_Win32Defs_H
