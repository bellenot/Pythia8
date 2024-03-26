#include "Pythia8/Win32Defs.h"

char *dlerror() {
   static char Msg[1000];
   FormatMessage(FORMAT_MESSAGE_FROM_SYSTEM, NULL, GetLastError(),
                 MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), Msg,
                 sizeof(Msg), NULL);
   return Msg;
}

FARPROC dlsym(void *library, const char *function_name)
{
   HMODULE hMods[1024];
   DWORD cbNeeded;
   FARPROC address = NULL;
   unsigned int i;
   if (library == RTLD_DEFAULT) {
      if (EnumProcessModules(::GetCurrentProcess(), hMods, sizeof(hMods), &cbNeeded)) {
         for (i = 0; i < (cbNeeded / sizeof(HMODULE)); i++) {
            address = ::GetProcAddress((HMODULE)hMods[i], function_name);
            if (address)
               return address;
         }
      }
      return address;
   } else {
      return ::GetProcAddress((HMODULE)library, function_name);
   }
}
