---
source: /Users/byron/projects/rockbeatspaper/src/routes/login/index.tsx
relative: rockbeatspaper/src/routes/login/index.tsx
generated_at: 2025-12-23 10:28
---

```tsx
import { createFileRoute, useNavigate, Link } from "@tanstack/react-router";
import { useForm } from "react-hook-form";
import { useMutation } from "@tanstack/react-query";
import toast from "react-hot-toast";
import { Navbar } from "~/components/Navbar";
import { Button } from "~/components/Button";
import { Input } from "~/components/Input";
import { Card } from "~/components/Card";
import { useTRPC } from "~/trpc/react";
import { useAuthStore } from "~/stores/auth";

export const Route = createFileRoute("/login/")({
  component: Login,
});

function Login() {
  const navigate = useNavigate();
  const trpc = useTRPC();
  const { setAuth } = useAuthStore();

  const {
    register,
    handleSubmit,
    formState: { errors },
  } = useForm<{ email: string; password: string }>();

  const loginMutation = useMutation(
    trpc.login.mutationOptions({
      onSuccess: (data) => {
        setAuth(data.authToken, data.user);
        toast.success("Welcome back!");
        navigate({ to: "/dashboard" });
      },
      onError: (error) => {
        toast.error(error.message || "Login failed");
      },
    })
  );

  const onSubmit = (data: { email: string; password: string }) => {
    loginMutation.mutate(data);
  };

  return (
    <div className="min-h-screen bg-gradient-to-br from-blue-50 via-purple-50 to-pink-50 dark:from-gray-900 dark:via-gray-800 dark:to-gray-900">
      <Navbar />
      <div className="flex items-center justify-center px-4 py-12">
        <div className="w-full max-w-md">
          <Card>
            <div className="text-center mb-8">
              <h1 className="text-3xl font-bold text-gray-900 dark:text-white mb-2">
                Welcome Back
              </h1>
              <p className="text-gray-600 dark:text-gray-400">
                Sign in to your account to continue
              </p>
            </div>

            <form onSubmit={handleSubmit(onSubmit)} className="space-y-6">
              <Input
                label="Email"
                type="email"
                {...register("email", { required: "Email is required" })}
                error={errors.email?.message}
              />

              <Input
                label="Password"
                type="password"
                {...register("password", { required: "Password is required" })}
                error={errors.password?.message}
              />

              <Button
                type="submit"
                disabled={loginMutation.isPending}
                className="w-full"
              >
                {loginMutation.isPending ? "Signing in..." : "Sign In"}
              </Button>
            </form>

            <div className="mt-6 text-center text-sm text-gray-600 dark:text-gray-400">
              Don't have an account?{" "}
              <Link
                to="/signup"
                className="text-blue-600 hover:text-blue-700 font-semibold"
              >
                Sign up
              </Link>
            </div>
          </Card>
        </div>
      </div>
    </div>
  );
}

```
