---
source: /Users/byron/projects/rockbeatspaper/src/routes/signup/index.tsx
relative: rockbeatspaper/src/routes/signup/index.tsx
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

export const Route = createFileRoute("/signup/")({
  component: Signup,
});

function Signup() {
  const navigate = useNavigate();
  const trpc = useTRPC();
  const { setAuth } = useAuthStore();

  const {
    register,
    handleSubmit,
    formState: { errors },
  } = useForm<{ name: string; email: string; password: string }>();

  const signupMutation = useMutation(
    trpc.signup.mutationOptions({
      onSuccess: (data) => {
        setAuth(data.authToken, data.user);
        toast.success("Account created successfully!");
        navigate({ to: "/dashboard" });
      },
      onError: (error) => {
        toast.error(error.message || "Signup failed");
      },
    })
  );

  const onSubmit = (data: { name: string; email: string; password: string }) => {
    signupMutation.mutate(data);
  };

  return (
    <div className="min-h-screen bg-gradient-to-br from-blue-50 via-purple-50 to-pink-50 dark:from-gray-900 dark:via-gray-800 dark:to-gray-900">
      <Navbar />
      <div className="flex items-center justify-center px-4 py-12">
        <div className="w-full max-w-md">
          <Card>
            <div className="text-center mb-8">
              <h1 className="text-3xl font-bold text-gray-900 dark:text-white mb-2">
                Create Account
              </h1>
              <p className="text-gray-600 dark:text-gray-400">
                Start optimizing your energy today
              </p>
            </div>

            <form onSubmit={handleSubmit(onSubmit)} className="space-y-6">
              <Input
                label="Full Name"
                type="text"
                {...register("name", { required: "Name is required" })}
                error={errors.name?.message}
              />

              <Input
                label="Email"
                type="email"
                {...register("email", { required: "Email is required" })}
                error={errors.email?.message}
              />

              <Input
                label="Password"
                type="password"
                {...register("password", {
                  required: "Password is required",
                  minLength: {
                    value: 8,
                    message: "Password must be at least 8 characters",
                  },
                })}
                error={errors.password?.message}
              />

              <Button
                type="submit"
                disabled={signupMutation.isPending}
                className="w-full"
              >
                {signupMutation.isPending ? "Creating account..." : "Sign Up"}
              </Button>
            </form>

            <div className="mt-6 text-center text-sm text-gray-600 dark:text-gray-400">
              Already have an account?{" "}
              <Link
                to="/login"
                className="text-blue-600 hover:text-blue-700 font-semibold"
              >
                Sign in
              </Link>
            </div>
          </Card>
        </div>
      </div>
    </div>
  );
}

```
