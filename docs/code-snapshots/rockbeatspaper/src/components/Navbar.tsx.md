---
source: /Users/byron/projects/rockbeatspaper/src/components/Navbar.tsx
relative: rockbeatspaper/src/components/Navbar.tsx
generated_at: 2025-12-23 10:28
---

```tsx
import { Link } from "@tanstack/react-router";
import { Zap, LogOut, User } from "lucide-react";
import { useAuthStore } from "~/stores/auth";
import { Button } from "./Button";

export function Navbar() {
  const { user, clearAuth } = useAuthStore();

  return (
    <nav className="bg-white dark:bg-gray-900 border-b border-gray-200 dark:border-gray-800 sticky top-0 z-50 backdrop-blur-lg bg-opacity-90">
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
        <div className="flex justify-between items-center h-16">
          <Link
            to="/"
            className="flex items-center space-x-2 text-2xl font-bold bg-gradient-to-r from-blue-600 to-purple-600 bg-clip-text text-transparent"
          >
            <Zap className="w-8 h-8 text-blue-600" />
            <span>ThermalAI</span>
          </Link>

          <div className="flex items-center space-x-4">
            {user ? (
              <>
                <Link to="/dashboard">
                  <Button variant="outline" size="sm">
                    Dashboard
                  </Button>
                </Link>
                <Link to="/pi-benford">
                  <Button variant="outline" size="sm">
                    Pi Analysis
                  </Button>
                </Link>
                <div className="flex items-center space-x-3 pl-4 border-l border-gray-300 dark:border-gray-700">
                  <div className="flex items-center space-x-2">
                    <div className="w-8 h-8 bg-gradient-to-r from-blue-600 to-purple-600 rounded-full flex items-center justify-center">
                      <User className="w-5 h-5 text-white" />
                    </div>
                    <span className="text-sm font-medium text-gray-700 dark:text-gray-300">
                      {user.name}
                    </span>
                  </div>
                  <button
                    onClick={clearAuth}
                    className="text-gray-600 hover:text-gray-900 dark:text-gray-400 dark:hover:text-gray-200"
                    title="Logout"
                  >
                    <LogOut className="w-5 h-5" />
                  </button>
                </div>
              </>
            ) : (
              <>
                <Link to="/pi-benford">
                  <Button variant="outline" size="sm">
                    Pi Analysis
                  </Button>
                </Link>
                <Link to="/login">
                  <Button variant="outline" size="sm">
                    Login
                  </Button>
                </Link>
                <Link to="/signup">
                  <Button variant="primary" size="sm">
                    Sign Up
                  </Button>
                </Link>
              </>
            )}
          </div>
        </div>
      </div>
    </nav>
  );
}

```
